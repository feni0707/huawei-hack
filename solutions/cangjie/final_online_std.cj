// final_online_std — порт solutions/final_online.cpp (= scorefit_guards.cpp: scorefit + предохранитель зон + страховка через
// тень best fit) на Cangjie (仓颉) 1.0. Стратегия и комментарии — в scorefit_guards.cpp / scorefit.cpp,
// здесь построчный перенос той же логики (имена функций и полей те же). Ввод-вывод: см. main().
//
// Базовый порт scorefit.cpp из коммита 4ff4b57 (до 3 кандидатов каждой длины в большой траншее). Проверено на
// cjc 1.0.5 (2026-10-20): компилируется без предупреждений, ответы побитово = C++ на materials/tests и build/gen/*;
// bench.py: open 30.1539 (= C++), 0 WA/TL.
// Сборка на macOS: `source build/cangjie-sdk/cangjie/envsetup.sh`, затем
//   SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk cjc -O2 scorefit.cj -o scorefit
// (lld из SDK не читает MacOSX26.sdk). Сравнить ответы с C++-версией:
//   for f in materials/tests/*.in build/gen/*/*.in; do tail -n +2 $f > /tmp/c.txt; ./scorefit < /tmp/c.txt > /tmp/a.txt;
//     build/bin/scorefit < /tmp/c.txt > /tmp/b.txt; cmp /tmp/a.txt /tmp/b.txt || echo DIFF $f; done
// Ответы должны совпадать побитово (логика и порядок операций с плавающей точкой — как в C++; std.math.log
// побитово равен libm log). Двоичный файл требует envsetup.sh и при запуске (DYLD_LIBRARY_PATH; из bash-скриптов
// macOS его вырезает — запускать через обёртку, которая сама делает source). Оговорка: clang на arm64 по умолчанию
// сливает a*b+c в FMA, cjc — нет; на ~4100 тестах из build/ это дало 2 расхождения, и оба исчезают при сравнении
// с C++, собранной с -ffp-contract=off (на x86-64 g++ без -march FMA нет, так что ответы как у этого порта).
//
// Отличия от C++ только технические: все целые — Int64 (страницы — Int32); массивы по PID — фиксированного
// размера MAXPID (команд ≤ 200 000); «нулевая» страница ленивого массива — пустой Array; тень — Option<Alloc>.
package final_online_std

import std.env.*
import std.math.*

let NB: Int64 = 128          // промежутки короче NB — в корзинах по точной длине, длиннее — в общей куче
let MAXL: Int64 = 40         // максимальная длина трубы
let CHUNK: Int64 = 400
let CHUNK_SAME: Int64 = 800
let CHUNK_MIN: Int64 = 200
let CHUNK_PER_L: Int64 = 10
let CHUNK_LOW: Int64 = 40
let CHUNK_LOW_AGE: Int64 = 1000
let EARLY_RUN: Int64 = 64
let EARLY_MAXL: Int64 = 3
let EARLY_CHUNK: Int64 = 40
let FIT_EXTRA: Int64 = 4
let FIT_RUN: Int64 = 30
let MULTI_K: Int64 = 3
let MULTI_T: Int64 = 50000
let SHADOW_MAX_T: Int64 = 3000000
let ZSLOTS: Int64 = 16
let ZSHARE: Float64 = 0.08
let ZMAXACT: Int64 = 8
let ZMAXSHARE: Float64 = 0.25
let ZMINALIVE: Int64 = 150
let ZFAST: Float64 = 0.7
let EXWARM: Float64 = 300.0
let SHIFT_TV: Float64 = 0.7
let SHIFT_GAP: Int64 = 500
let MSEEN: Int64 = 1000
let ZLIFO: Float64 = 0.7
let GUARD_PCT: Float64 = 10.0      // страховка: наш T > T_тени·(1 + GUARD_PCT/100) и T − T_тени ≥ GUARD_ABS → best fit
let GUARD_ABS: Int64 = 1000
let GUARD_MIN: Int64 = 2000
let ZGUARD_FREE: Float64 = 0.15    // предохранитель зон: доля свободного, выше которой зоны не копают
let NF: Int64 = 14
let BITN: Int64 = 262144     // 1 << 18 > 200 000 укладок
let MAXPID: Int64 = 200002
let PSH: Int64 = 16          // страница 2^16
let PMASK: Int64 = 65535

// веса линейной оценки кандидата (см. scorefit.cpp)
let LW: Array<Float64> = [1.0, -38.8, -22.4, 13.5, 47.0, -203.0, 107.0, -7.25, 2.78, 2.26, -8.96, 23.4, 1050.0, 69.0,
    0.0, 46.7, 27.9, 50.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 50.0]

// таблицы логарифмов (те же значения, что log, без вызова в горячем цикле)
let LOGN: Int64 = 4096
let logN: Array<Float64> = Array<Float64>(LOGN, {n => log(1.0 + Float64(n))})
let logSL: Array<Float64> = Array<Float64>(NB * (MAXL + 1), {i =>
    let s = i / (MAXL + 1)
    let l = i % (MAXL + 1)
    if (s > 0 && l >= 1) { log(Float64(s) / Float64(l)) } else { 0.0 }
})
let EMPTY: Array<Int32> = Array<Int32>()

func dmax(a: Float64, b: Float64): Float64 { if (a > b) { a } else { b } }
func dmin(a: Float64, b: Float64): Float64 { if (a < b) { a } else { b } }
func dabs(a: Float64): Float64 { if (a < 0.0) { -a } else { a } }
func imax(a: Int64, b: Int64): Int64 { if (a > b) { a } else { b } }
func b2f(b: Bool): Float64 { if (b) { 1.0 } else { 0.0 } }
func log1n(n: Int64): Float64 { if (n < LOGN) { logN[n] } else { log(1.0 + Float64(n)) } }
// верхняя граница разброса gapCost между промежутками одной длины (признаки соседей: LW7, 8, 9, 15, 17)
func nbW(r: Int64): Float64 {
    dabs(LW[7]) + dabs(LW[8]) * Float64(r) + dabs(LW[9]) + dabs(LW[15]) + dabs(LW[17]) + 1.0
}
func ageBucket(a: Int64): Int64 {
    var b: Int64 = 0
    while ((2 << b) - 1 <= a && b < 23) { b += 1 }
    return b
}

// ---------------------------------------------------------------- растущий массив Int64 (стек)
class IVec {
    var a: Array<Int64>
    var n: Int64
    init() { a = Array<Int64>(4, repeat: 0); n = 0 }
    init(o: IVec) { a = o.a.clone(); n = o.n }
    func size(): Int64 { n }
    func empty(): Bool { n == 0 }
    func get(i: Int64): Int64 { a[i] }
    func set(i: Int64, v: Int64): Unit { a[i] = v }
    func back(): Int64 { a[n - 1] }
    func clear(): Unit { n = 0 }
    func push(v: Int64): Unit {
        if (n == a.size) {
            let b = Array<Int64>(a.size * 2, repeat: 0)
            a.copyTo(b, 0, 0, n)
            a = b
        }
        a[n] = v
        n += 1
    }
    func pop(): Int64 { n -= 1; a[n] }
}

// ---------------------------------------------------------------- постраничный массив (Int32)
// Ленивый режим: страница создаётся при первой ненулевой записи и освобождается, когда снова вся нулевая.
class Paged {
    var pg: Array<Array<Int32>>
    var nz: Array<Int64>
    var npg: Int64
    let lazy: Bool
    let spare: Array<Array<Int32>> = Array<Array<Int32>>(4, repeat: EMPTY)
    var nsp: Int64 = 0
    init(lz: Bool) { lazy = lz; pg = Array<Array<Int32>>(16, repeat: EMPTY); nz = Array<Int64>(16, repeat: 0); npg = 0 }
    init(o: Paged) {   // глубокая копия (тень аллокатора)
        lazy = o.lazy
        npg = o.npg
        nz = o.nz.clone()
        let op = o.pg
        let on = o.npg
        pg = Array<Array<Int32>>(op.size, {k => if (k < on && op[k].size > 0) { op[k].clone() } else { EMPTY }})
    }
    func newPage(): Array<Int32> {
        if (nsp > 0) { nsp -= 1; let p = spare[nsp]; spare[nsp] = EMPTY; return p }
        return Array<Int32>(PMASK + 1, repeat: 0)
    }
    func grow(n: Int64): Unit {
        while ((npg << PSH) < n) {
            if (npg == pg.size) {
                let b = Array<Array<Int32>>(pg.size * 2, repeat: EMPTY)
                pg.copyTo(b, 0, 0, npg)
                pg = b
                let c = Array<Int64>(nz.size * 2, repeat: 0)
                nz.copyTo(c, 0, 0, npg)
                nz = c
            }
            pg[npg] = if (lazy) { EMPTY } else { newPage() }
            nz[npg] = 0
            npg += 1
        }
    }
    func get(i: Int64): Int64 {
        let p = pg[i >> PSH]
        if (p.size == 0) { 0 } else { Int64(p[i & PMASK]) }
    }
    func set(i: Int64, v: Int64): Unit {
        let k = i >> PSH
        if (!lazy) { pg[k][i & PMASK] = Int32(v); return }
        var p = pg[k]
        if (p.size == 0) {
            if (v == 0) { return }
            p = newPage()
            pg[k] = p
        }
        let c = Int64(p[i & PMASK])
        nz[k] += (if (v != 0) { 1 } else { 0 }) - (if (c != 0) { 1 } else { 0 })
        p[i & PMASK] = Int32(v)
        if (nz[k] == 0) {
            if (nsp < 4) { spare[nsp] = p; nsp += 1 }
            pg[k] = EMPTY
        }
    }
}

// ---------------------------------------------------------------- двоичная куча (минимум)
class Heap {
    var a: Array<Int64>
    var n: Int64
    init() { a = Array<Int64>(16, repeat: 0); n = 0 }
    init(o: Heap) { a = o.a.clone(); n = o.n }
    func empty(): Bool { n == 0 }
    func top(): Int64 { a[0] }
    func push(v: Int64): Unit {
        if (n == a.size) {
            let b = Array<Int64>(a.size * 2, repeat: 0)
            a.copyTo(b, 0, 0, n)
            a = b
        }
        a[n] = v
        var i = n
        n += 1
        while (i > 0) {
            let p = (i - 1) / 2
            if (a[p] <= a[i]) { break }
            let t = a[p]; a[p] = a[i]; a[i] = t
            i = p
        }
    }
    func pop(): Unit {
        n -= 1
        a[0] = a[n]
        var i: Int64 = 0
        while (true) {
            let l = 2 * i + 1
            let r = l + 1
            var m = i
            if (l < n && a[l] < a[m]) { m = l }
            if (r < n && a[r] < a[m]) { m = r }
            if (m == i) { break }
            let t = a[m]; a[m] = a[i]; a[i] = t
            i = m
        }
    }
}

// ---------------------------------------------------------------- аллокатор
class Alloc {
    var T: Int64 = 0
    var tag: Paged
    var owner: Paged
    var zl: Paged
    var bucket: Array<Heap>
    var big: Heap
    var gapCnt: Array<Int64>

    var npid: Int64 = 1                  // = pos.size() в C++
    var pos: Array<Int64>
    var len: Array<Int64>
    var alive: Array<Bool>
    var feat: Array<Int64>
    var stream: Array<Int64>

    var cnt: Array<Float64>
    var cntTot: Float64 = 0.0
    var lastSeen: Array<Int64>
    var cntR: Array<Float64>
    var lastShift: Int64 = 0
    var birthL: Array<Float64>
    var deathL: Array<Float64>
    var recDel: Array<Float64>
    var recDelTot: Float64 = 0.0
    var aliveCntL: Array<Int64>
    var aliveTot: Int64 = 0

    var F: Int64 = 0
    var flip: Bool = false
    var noStreams: Bool = false
    var shadowDropped: Bool = false
    var shadow: Option<Alloc> = None
    var shadowBorn: Int64 = 0
    var pureBF: Bool = false             // сама является тенью чистого best fit (без потоков и зон)
    var gsh: Option<Alloc> = None        // тень чистого best fit на той же истории команд
    var gDropped: Bool = false
    var gOn: Bool = false                // основной сейчас решает как best fit
    var occ: Int64 = 0
    var runIdx: Int64 = 0
    var lastL: Int64 = 0
    var runPrev: Int64 = 0
    var runLast: Int64 = 0
    var coSame: Array<Float64>
    var coTot: Array<Float64>
    var birthV: Array<Float64>           // [f * 2 + v]
    var deathV: Array<Float64>
    var curFeat: Int64 = -1
    var bornLab: Array<Float64>
    var diedLab: Array<Float64>
    var qSum: Float64 = 0.0
    var qTot: Float64 = 0.0
    var raceOld: Float64 = 0.0
    var raceTot: Float64 = 0.0
    var raceHead: Array<Int64>
    var raceNext: IVec
    var raceOther: IVec
    var raceDone: IVec
    var lastDelPid: Int64 = 0
    var lastWasDel: Bool = false

    var aliveBit: Array<Int64>
    var deathAge: Array<Float64>
    var deathSuf: Array<Float64>
    var anyDel: Bool = false
    var verB: Bool = false

    var slots: Array<IVec>
    var resEnd: Int64 = 0
    var cntF: Array<Float64>
    var cntFTot: Float64 = 0.0
    var lastOfLen: Array<Int64>          // [l * 8 + k]
    var lastPtr: Array<Int64>
    var lifoHit: Float64 = 0.0
    var lifoTot: Float64 = 0.0

    // временные (не копируются)
    let tmpKeep: IVec = IVec()
    var outLen: Int64 = 0                // длина промежутка из firstBig / scanBig
    var zN: Int64 = 0                    // результат zoneTry
    var zD: Int64 = -1
    let gk: Array<Int64> = Array<Int64>(MULTI_K, repeat: 0)

    init() {
        tag = Paged(true); owner = Paged(false); zl = Paged(true)
        bucket = Array<Heap>(NB, {_ => Heap()})
        big = Heap()
        gapCnt = Array<Int64>(NB, repeat: 0)
        pos = Array<Int64>(MAXPID, repeat: 0)
        pos[0] = -1
        len = Array<Int64>(MAXPID, repeat: 0)
        alive = Array<Bool>(MAXPID, repeat: false)
        feat = Array<Int64>(MAXPID, repeat: 0)
        stream = Array<Int64>(MAXPID, repeat: -1)
        cnt = Array<Float64>(MAXL + 1, repeat: 0.0)
        lastSeen = Array<Int64>(MAXL + 1, repeat: 0)
        cntR = Array<Float64>(MAXL + 1, repeat: 0.0)
        birthL = Array<Float64>(MAXL + 1, repeat: 0.0)
        deathL = Array<Float64>(MAXL + 1, repeat: 0.0)
        recDel = Array<Float64>(MAXL + 1, repeat: 0.0)
        aliveCntL = Array<Int64>(MAXL + 1, repeat: 0)
        coSame = Array<Float64>(NF, repeat: 0.0)
        coTot = Array<Float64>(NF, repeat: 0.0)
        birthV = Array<Float64>(NF * 2, repeat: 0.0)
        deathV = Array<Float64>(NF * 2, repeat: 0.0)
        bornLab = Array<Float64>(2, repeat: 0.0)
        diedLab = Array<Float64>(2, repeat: 0.0)
        raceHead = Array<Int64>(MAXPID, repeat: -1)
        raceNext = IVec(); raceOther = IVec(); raceDone = IVec()
        aliveBit = Array<Int64>(BITN, repeat: 0)
        deathAge = Array<Float64>(24, repeat: 0.0)
        deathSuf = Array<Float64>(24, repeat: 0.0)
        slots = Array<IVec>(MAXL + 1, {_ => IVec()})
        cntF = Array<Float64>(MAXL + 1, repeat: 0.0)
        lastOfLen = Array<Int64>((MAXL + 1) * 8, repeat: 0)
        lastPtr = Array<Int64>(MAXL + 1, repeat: 0)
        ensure(1024)
    }
    init(o: Alloc) {   // копия для тени (shadow не копируется)
        T = o.T
        tag = Paged(o.tag); owner = Paged(o.owner); zl = Paged(o.zl)
        let ob = o.bucket
        bucket = Array<Heap>(NB, {i => Heap(ob[i])})
        big = Heap(o.big)
        gapCnt = o.gapCnt.clone()
        npid = o.npid
        pos = o.pos.clone(); len = o.len.clone(); alive = o.alive.clone(); feat = o.feat.clone(); stream = o.stream.clone()
        cnt = o.cnt.clone(); cntTot = o.cntTot; lastSeen = o.lastSeen.clone(); cntR = o.cntR.clone(); lastShift = o.lastShift
        birthL = o.birthL.clone(); deathL = o.deathL.clone(); recDel = o.recDel.clone(); recDelTot = o.recDelTot
        aliveCntL = o.aliveCntL.clone(); aliveTot = o.aliveTot
        F = o.F; flip = o.flip; noStreams = o.noStreams; shadowDropped = o.shadowDropped; shadowBorn = o.shadowBorn
        pureBF = o.pureBF; gDropped = o.gDropped; gOn = o.gOn
        occ = o.occ; runIdx = o.runIdx; lastL = o.lastL; runPrev = o.runPrev; runLast = o.runLast
        coSame = o.coSame.clone(); coTot = o.coTot.clone(); birthV = o.birthV.clone(); deathV = o.deathV.clone()
        curFeat = o.curFeat; bornLab = o.bornLab.clone(); diedLab = o.diedLab.clone()
        qSum = o.qSum; qTot = o.qTot; raceOld = o.raceOld; raceTot = o.raceTot
        raceHead = o.raceHead.clone(); raceNext = IVec(o.raceNext); raceOther = IVec(o.raceOther); raceDone = IVec(o.raceDone)
        lastDelPid = o.lastDelPid; lastWasDel = o.lastWasDel
        aliveBit = o.aliveBit.clone(); deathAge = o.deathAge.clone(); deathSuf = o.deathSuf.clone()
        anyDel = o.anyDel; verB = o.verB
        let os = o.slots
        slots = Array<IVec>(MAXL + 1, {i => IVec(os[i])})
        resEnd = o.resEnd; cntF = o.cntF.clone(); cntFTot = o.cntFTot
        lastOfLen = o.lastOfLen.clone(); lastPtr = o.lastPtr.clone(); lifoHit = o.lifoHit; lifoTot = o.lifoTot
    }

    // ---- промежутки
    func ensure(n: Int64): Unit { tag.grow(n + 1); owner.grow(n + 1); zl.grow(n + 1) }
    func addGap(s: Int64, l: Int64): Unit {
        tag.set(s, l); tag.set(s + l, -s - 1)
        if (l < NB) { bucket[l].push(s); gapCnt[l] += 1 } else { big.push((l << 32) | s) }
    }
    func delGap(s: Int64): Unit {
        let l = tag.get(s)
        if (l < NB) { gapCnt[l] -= 1 }
        tag.set(s, 0); tag.set(s + l, 0)
    }
    func freeRange(s0: Int64, l: Int64): Int64 {   // освободить [s, s+l) со склейкой соседей; вернуть конец
        var s = s0
        var e = s0 + l
        if (e < T && tag.get(e) > 0) { let r = tag.get(e); delGap(e); e += r }
        if (tag.get(s) < 0) { let ps = -tag.get(s) - 1; delGap(ps); s = ps }
        addGap(s, e - s)
        return e
    }
    func takeFrom(gs: Int64, d: Int64, l: Int64): Unit {
        let gl = tag.get(gs)
        delGap(gs)
        if (d > gs) { addGap(gs, d - gs) }
        if (d + l < gs + gl) { addGap(d + l, gs + gl - (d + l)) }
    }
    func dig(n: Int64): Unit {
        if (n <= 0) { return }
        let oldT = T
        T += 20 * n
        ensure(T)
        freeRange(oldT, T - oldT)
    }
    func tailStart(): Int64 { if (T > 0 && tag.get(T) < 0) { -tag.get(T) - 1 } else { T } }
    func pidLeftOf(x: Int64): Int64 { if (x > 0) { owner.get(x - 1) } else { 0 } }
    func pidRightOf(x: Int64): Int64 { if (x < T) { owner.get(x) } else { 0 } }

    func firstGap(s: Int64, ts: Int64): Int64 {
        if (gapCnt[s] == 0) { return -1 }
        let h = bucket[s]
        var k0: Int64 = 0
        var k1: Int64 = 0
        var k2: Int64 = 0
        var nk: Int64 = 0
        var res: Int64 = -1
        while (!h.empty() && nk < 3) {
            let g = h.top()
            let lastK = if (nk == 1) { k0 } else if (nk == 2) { k1 } else { -1 }
            if (tag.get(g) != s || (nk > 0 && lastK == g)) { h.pop(); continue }
            if (g != ts && g + s != F) { res = g; break }
            h.pop()
            if (nk == 0) { k0 = g } else if (nk == 1) { k1 = g } else { k2 = g }
            nk += 1
        }
        if (nk > 0) { h.push(k0) }
        if (nk > 1) { h.push(k1) }
        if (nk > 2) { h.push(k2) }
        return res
    }
    // до k действительных промежутков длины s (по адресу), кроме хвоста и граничного; элементы кучи — на место
    func someGaps(s: Int64, ts: Int64, k: Int64, out: Array<Int64>): Int64 {
        if (gapCnt[s] == 0) { return 0 }
        let h = bucket[s]
        tmpKeep.clear()
        var no: Int64 = 0
        while (!h.empty() && tmpKeep.size() < 16 && no < k) {
            let g = h.top()
            h.pop()
            if (tag.get(g) != s || (tmpKeep.size() > 0 && tmpKeep.back() == g)) { continue }
            tmpKeep.push(g)
            if (g != ts && g + s != F) { out[no] = g; no += 1 }
        }
        for (i in 0..tmpKeep.size()) { h.push(tmpKeep.get(i)) }
        return no
    }
    // самый короткий большой промежуток (среди равных — левый); длина — в outLen
    func firstBig(ts: Int64): Int64 {
        tmpKeep.clear()
        var res: Int64 = -1
        while (!big.empty() && tmpKeep.size() < 3) {
            let k = big.top()
            let l = k >> 32
            let g = k & 0xFFFFFFFF
            if (tag.get(g) != l || (tmpKeep.size() > 0 && tmpKeep.back() == k)) { big.pop(); continue }
            big.pop()
            tmpKeep.push(k)
            if (g != ts && g + l != F) { res = g; outLen = l; break }
        }
        for (i in 0..tmpKeep.size()) { big.push(tmpKeep.get(i)) }
        return res
    }

    // ---- статистика
    func minRecent(): Int64 {
        var a: Int64 = 1
        for (l in 1..=MAXL) { if (cnt[l] > 0.002 * cntTot) { a = l; break } }
        let np = npid - 1
        for (l in 1..=MAXL) {
            if (lastSeen[l] > 0 && np - lastSeen[l] < MSEEN) { if (l > a) { a = l }; break }
        }
        return a
    }
    func pdieLen(p: Int64): Float64 {
        if (p <= 0) { return 0.0 }
        let l = len[p]
        return (deathL[l] + 0.5) / (birthL[l] + 1.0)
    }
    func rateOfLen(l: Int64): Float64 {
        let rel = (recDel[l] / recDelTot) / (Float64(imax(1, aliveCntL[l])) / Float64(imax(1, aliveTot)))
        return dmin(1.0, rel / 4.0)
    }
    func hzRate(p: Int64): Float64 { if (p <= 0 || recDelTot < 1.0) { 0.0 } else { rateOfLen(len[p]) } }
    func demand(z: Int64): Float64 {
        if (z > MAXL) { 1e9 } else if (cntTot > 0.0) { cnt[z] / cntTot * LW[26] } else { 0.0 }
    }

    func bitAdd(i0: Int64, v: Int64): Unit {
        var i = i0
        while (i < BITN) { aliveBit[i] += v; i += i & -i }
    }
    func bitSum(i0: Int64): Int64 {
        var i = i0
        var r: Int64 = 0
        while (i > 0) { r += aliveBit[i]; i -= i & -i }
        return r
    }
    func hazard(p: Int64): Float64 {
        if (p <= 0) { return 0.0 }
        let cur = npid - 1
        let b = ageBucket(cur - p)
        let dGE = deathSuf[b]
        let lo = cur - ((1 << b) - 1)
        let aGE = if (lo > 0) { Float64(bitSum(lo)) } else { 0.0 }
        return if (dGE + aGE > 0.0) { dGE / (dGE + aGE) } else { 0.5 }
    }

    // ---- потоки
    func kappa(f: Int64): Float64 {
        if (coTot[f] < 20.0) { return 0.0 }
        let b0 = birthV[f * 2]
        let b1 = birthV[f * 2 + 1]
        let bt = b0 + b1
        if (bt <= 0.0) { return 0.0 }
        if (f >= 5 && runLast < FIT_RUN) { return 0.0 }
        let minSh = if (f < 5) { 0.3 } else if (f < 9) { 0.2 } else { 0.16 }
        let mn = if (b0 < b1) { b0 } else { b1 }
        if (mn < minSh * bt) { return 0.0 }
        let e = (b0 * b0 + b1 * b1) / (bt * bt)
        let o = coSame[f] / coTot[f]
        return if (e >= 1.0) { 0.0 } else { (o - e) / (1.0 - e) }
    }
    func pdie(f: Int64, v: Int64): Float64 { (deathV[f * 2 + v] + 1.0) / (birthV[f * 2 + v] + 2.0) }
    func chooseStream(pid: Int64): Int64 {
        if (noStreams) { return -1 }
        if (!anyDel && runIdx >= EARLY_RUN && lastL <= EARLY_MAXL) { return runIdx & 1 }
        var bf: Int64 = -1
        var bs: Float64 = -1.0
        var cs: Float64 = -1.0
        for (f in 0..NF) {
            let k = kappa(f)
            if (k < 0.6) { continue }
            let s = k + dabs(pdie(f, 0) - pdie(f, 1))
            if (f == curFeat) { cs = s }
            if (s > bs) { bs = s; bf = f }
        }
        if (bf < 0) { curFeat = -1; return -1 }
        if (cs >= 0.0 && bs < cs + 0.2) { bf = curFeat }
        curFeat = bf
        let v = (feat[pid] >> bf) & 1
        let pv = pdie(bf, v)
        let po = pdie(bf, v ^ 1)
        if (dabs(pv - po) < 0.1) { return v }
        return if (pv > po) { 1 } else { 0 }
    }

    func zoneBlocked(gs: Int64, s: Int64, st: Int64): Bool {
        if (!verB || st < 0) { return false }
        let lp = pidLeftOf(gs)
        let rp = pidRightOf(gs + s)
        if (st == 1) { return (lp > 0 && stream[lp] != 1) || (rp > 0 && stream[rp] != 1) }
        return lp > 0 && rp > 0 && stream[lp] == 1 && stream[rp] == 1
    }

    func chooseSide(gs: Int64, s: Int64, L: Int64, st: Int64): Int64 {
        let r = s - L
        if (r == 0) { return gs }
        let lp = pidLeftOf(gs)
        let rp = pidRightOf(gs + s)
        if (st >= 0) {
            let lc = if (lp != 0) { stream[lp] } else { -2 }
            let rc = if (rp != 0) { stream[rp] } else { -2 }
            if (lc == st && rc != st) { return gs }
            if (rc == st && lc != st) { return gs + r }
            return if (st == 1) { gs + r } else { gs }
        }
        return if (sideNb(lp, rp) == lp) { gs + r } else { gs }
    }
    func oldDieFirst(): Bool {
        if (raceTot >= 20.0 && dabs(raceOld / raceTot - 0.5) > 0.2) { return raceOld > 0.5 * raceTot }
        return qTot >= 20.0 && qSum / qTot > 0.5
    }
    func raceAdd(a: Int64, b: Int64): Unit {
        let id = raceDone.size()
        raceDone.push(0)
        raceNext.push(raceHead[a]); raceOther.push(b); raceHead[a] = 2 * id
        raceNext.push(raceHead[b]); raceOther.push(a); raceHead[b] = 2 * id + 1
    }
    func raceOnDeath(p: Int64): Unit {
        var e = raceHead[p]
        while (e >= 0) {
            if (raceDone.get(e >> 1) == 0) {
                raceDone.set(e >> 1, 1)
                raceOld = raceOld * 0.99 + (if (p < raceOther.get(e)) { 1.0 } else { 0.0 })
                raceTot = raceTot * 0.99 + 1.0
            }
            e = raceNext.get(e)
        }
        raceHead[p] = -1
    }
    func sideNb(lp: Int64, rp: Int64): Int64 {
        if (oldDieFirst()) {
            if (lp == 0) { return rp }
            if (rp == 0) { return lp }
            return if (lp < rp) { lp } else { rp }
        }
        return if (lp > rp) { lp } else { rp }
    }

    // стоимость укладки трубы L в промежуток gs длины s (меньше — лучше)
    func gapCost(gs: Int64, s: Int64, L: Int64, st: Int64, m: Int64, selfHz: Float64): Float64 {
        if (T < 60) {
            if (s == L) { return -1000.0 + 2000.0 * dmax(pdieLen(pidLeftOf(gs)), pdieLen(pidRightOf(gs + s))) }
            return -Float64(s)
        }
        let r = s - L
        let lp = pidLeftOf(gs)
        let rp = pidRightOf(gs + s)
        var seal: Float64 = 0.0
        if (r == 0) {
            var h: Float64 = 0.0
            for (k in 0..2) {
                let nb = if (k == 0) { lp } else { rp }
                if (nb <= 0) { continue }
                let mg = s + len[nb]
                if (mg > MAXL || (cnt[mg] < 0.05 * cntTot && cntF[mg] < 0.2 * cntFTot)) { continue }
                h = dmax(h, hzRate(nb))
            }
            h -= selfHz
            if (h >= 0.5) { seal = h }
        }
        if (st >= 0) { return if (r != 0) { Float64(r) } else { 3.0 * seal } }
        let hl = hazard(lp)
        let hr = hazard(rp)
        let hside = if (sideNb(lp, rp) == lp) { hl } else { hr }
        let dead = r > 0 && r < m
        let dr = if (r >= 1 && r <= MAXL && cntTot > 0.0) { cntR[r] / cntTot } else { 0.0 }
        let ds = if (s <= MAXL && cntTot > 0.0) { cntR[s] / cntTot } else { 0.0 }
        let ns = if (s < NB) { Float64(gapCnt[s]) } else { 0.0 }
        let nr = if (r < NB) { Float64(gapCnt[r]) } else { 0.0 }
        let fr = Float64(r)
        let lsl = if (s < NB) { logSL[s * (MAXL + 1) + L] } else { log(Float64(s) / Float64(L)) }
        var c = LW[0] * fr + LW[1] * b2f(r == 0) + LW[2] * b2f(dead) + LW[3] * (if (dead) { fr } else { 0.0 }) +
            LW[4] * b2f(r >= m) + LW[5] * dr + LW[6] * ds + LW[7] * (if (r > 0) { hside } else { 0.0 }) +
            LW[8] * (if (dead) { fr * (1.0 - dmax(hl, hr)) } else { 0.0 }) +
            LW[9] * (if (r == 0) { dmax(hzRate(lp), hzRate(rp)) } else { 0.0 }) + LW[10] * log1n(Int64(ns)) +
            LW[11] * lsl + LW[12] * b2f(dead && r >= 4) + LW[15] * (if (r == 0) { dmin(hl, hr) } else { 0.0 }) +
            LW[16] * b2f(r > 0 && r < 2 * m && !dead)
        if (r > 0) { c += LW[18] * ds / (1.0 + ns) + LW[19] * dr / (1.0 + nr) }
        if (r > 0 && cntTot >= EXWARM) {
            let exr = nr + 1.0 - demand(r)
            let exs = ns - demand(s)
            c += LW[25] * (log(1.0 + dmax(0.0, exr)) - log(1.0 + dmax(0.0, exs)))
        }
        return c + LW[17] * seal
    }

    // ---- зоны по длинам
    func zonesEnabled(): Bool {
        if (cntTot < 50.0 || T < 60) { return false }
        var k: Int64 = 0
        var mx: Float64 = 0.0
        for (l in 1..=MAXL) {
            let sh = cnt[l] / cntTot
            if (sh > ZSHARE) { k += 1 }
            mx = dmax(mx, sh)
        }
        if (!(k >= 4 && k <= ZMAXACT && mx <= ZMAXSHARE)) { return false }
        if (lifoTot < 20.0 || lifoHit < ZLIFO * lifoTot) { return false }
        var f: Float64 = 0.0
        for (l in 1..=MAXL) { if (cnt[l] > ZSHARE * cntTot) { f += cntF[l] } }
        return f >= ZFAST * cntFTot
    }
    func zoneOn(z: Int64, en: Bool): Bool { en && z > 0 && cnt[z] > ZSHARE * cntTot }
    func occupied(x: Int64): Bool {
        let p = owner.get(x)
        return p > 0 && alive[p] && pos[p] <= x && x < pos[p] + len[p]
    }
    func clearOf(d: Int64, n: Int64, L: Int64, en: Bool): Bool {
        for (x in d..d + n) {
            let z = zl.get(x)
            if (z != 0 && z != L && zoneOn(z, en)) { return false }
        }
        return true
    }
    // mode: 0 — точное попадание в свою/незарезервированную зону, 1/2 — «редкая длина» (не хвост, не чужая зона)
    func scanGap(s: Int64, L: Int64, mode: Int64, ts: Int64, en: Bool): Int64 {
        let h = bucket[s]
        tmpKeep.clear()
        var res: Int64 = -1
        while (!h.empty() && tmpKeep.size() < 64) {
            let g = h.top()
            h.pop()
            if (tag.get(g) != s || (!tmpKeep.empty() && tmpKeep.back() == g)) { continue }
            tmpKeep.push(g)
            let ok = if (mode == 0) {
                zl.get(g) == L || !zoneOn(zl.get(g), en)
            } else {
                g != ts && !zoneOn(zl.get(g), en) && clearOf(g, if (mode == 2) { s } else { L }, L, en)
            }
            if (ok) { res = g; break }
        }
        for (i in 0..tmpKeep.size()) { h.push(tmpKeep.get(i)) }
        return res
    }
    func scanBig(L: Int64, ts: Int64, en: Bool, need: Int64): Int64 {   // длина — в outLen
        tmpKeep.clear()
        var res: Int64 = -1
        while (!big.empty() && tmpKeep.size() < 64) {
            let k = big.top()
            big.pop()
            let l = k >> 32
            let g = k & 0xFFFFFFFF
            if (tag.get(g) != l || (!tmpKeep.empty() && tmpKeep.back() == k)) { continue }
            tmpKeep.push(k)
            if (g != ts && !zoneOn(zl.get(g), en) && clearOf(g, if (need < l) { need } else { l }, L, en)) {
                res = g; outLen = l; break
            }
        }
        for (i in 0..tmpKeep.size()) { big.push(tmpKeep.get(i)) }
        return res
    }
    // остаток >= m самого редкого размера среди незарезервированных дыр (не хвост); -1 — нет такой
    func rarePlace(L: Int64, m: Int64, ts: Int64, en: Bool): Int64 {
        var bc: Float64 = 1e18
        var bg: Int64 = -1
        for (s in L + m..NB) {
            let g = scanGap(s, L, 1, ts, en)
            if (g < 0) { continue }
            let c = Float64(gapCnt[s - L]) * 1000.0 + Float64(s)
            if (c < bc) { bc = c; bg = g }
        }
        let g = scanBig(L, ts, en, L)
        if (g >= 0) {
            let c = Float64(outLen)
            if (c < bc) { bc = c; bg = g }
        }
        if (bg >= 0) { takeFrom(bg, bg, L); return bg }
        return -1
    }
    func inflated(): Bool { Float64(T - occ) > ZGUARD_FREE * Float64(T) }
    func zoneTry(L: Int64, st: Int64, m: Int64): Bool {   // результат — zN, zD
        if (st >= 0) { return false }
        let en = zonesEnabled()
        if (!en) { return false }
        let ts = tailStart()
        zN = 0
        zD = -1
        if (L < NB) {
            let g = scanGap(L, L, 0, ts, en)
            if (g >= 0) { zD = g; takeFrom(g, g, L); return true }
        }
        if (zoneOn(L, en) && aliveCntL[L] >= ZMINALIVE) {
            let sl = slots[L]
            while (!sl.empty()) {
                let d = sl.pop()
                if (zl.get(d) != L || occupied(d) || occupied(d + L - 1)) { continue }
                var hs = d
                var steps: Int64 = 0
                while (hs > 0 && !occupied(hs - 1) && steps < 4096) { hs -= 1; steps += 1 }
                if (steps >= 4096 || tag.get(hs) <= 0 || hs + tag.get(hs) < d + L) { continue }
                zD = d
                takeFrom(hs, d, L)
                return true
            }
            var hl: Int64 = 0
            var hg: Int64 = -1
            var s2 = 4 * L
            while (s2 < NB && hg < 0) { hg = scanGap(s2, L, 2, ts, en); s2 += 1 }
            if (hg >= 0) {
                hl = tag.get(hg)
            } else {
                hg = scanBig(L, ts, en, ZSLOTS * L)
                if (hg >= 0) { hl = outLen }
            }
            if (hg >= 0) {
                var C2 = ZSLOTS * L
                let fit = (hl / L) * L
                if (fit < C2) { C2 = fit }
                for (x in hg..hg + C2) { zl.set(x, L) }
                var x = hg + C2 - L
                while (x > hg) { sl.push(x); x -= L }
                zD = hg
                takeFrom(hg, hg, L)
                return true
            }
            let rp = rarePlace(L, m, ts, en)
            if (rp >= 0) { zD = rp; return true }
            let C = ZSLOTS * L
            let stp = if (ts > resEnd) { ts } else { resEnd }
            if (T < stp + C && inflated()) { return false }   // предохранитель: раздутая траншея — не копать под зону
            if (T < stp + C) { zN = (stp + C - T + 19) / 20; dig(zN) }
            for (x in stp..stp + C) { zl.set(x, L) }
            var x = stp + C - L
            while (x > stp) { sl.push(x); x -= L }
            resEnd = stp + C
            zD = stp
            takeFrom(tailStart(), stp, L)
            return true
        }
        let rp = rarePlace(L, m, ts, en)
        if (rp >= 0) { zD = rp; return true }
        let stp = if (ts > resEnd) { ts } else { resEnd }
        if (T < stp + L && inflated()) { return false }   // предохранитель: раздутая траншея — не копать под редкую длину
        if (T < stp + L) { zN = (stp + L - T + 19) / 20; dig(zN) }
        zD = stp
        takeFrom(tailStart(), stp, L)
        if (resEnd < stp + L && resEnd > ts) { resEnd = stp + L }
        return true
    }

    // ---- тень без потоков
    func freeIn(): Int64 {
        var f = T - occ - (T - tailStart())
        if (F < T && F > 0 && tag.get(F) < 0) { f -= F - (-tag.get(F) - 1) }
        return f
    }
    func shadowCheck(): Unit {
        if (shadowDropped || noStreams) { return }
        if (shadow.isSome()) {
            if (T > SHADOW_MAX_T) { shadow = None; shadowDropped = true }
            return
        }
        var any = false
        for (f in 0..NF) { if (kappa(f) >= 0.6) { any = true } }
        if (!any) { return }
        shadowBorn = npid
        let sh = Alloc(this)
        sh.noStreams = true
        shadow = Some(sh)
    }

    func fitChunk(tl: Int64, N: Int64, L: Int64): Int64 {
        if (runIdx < FIT_RUN) { return N }
        var best = N
        var br = (tl + 20 * N) % L
        var n = N + 1
        while (n <= N + FIT_EXTRA && br > 0) {
            let r = (tl + 20 * n) % L
            if (r < br) { br = r; best = n }
            n += 1
        }
        return best
    }

    // страховка: решать как основной, если его T больше тени чистого best fit
    func guardOn(): Bool {
        if (gDropped) { return false }
        let g = match (gsh) { case Some(x) => x; case None => return false }
        if (T > SHADOW_MAX_T) { gsh = None; gDropped = true; gOn = false; return false }
        if (npid - 1 < GUARD_MIN) { return false }
        let x = GUARD_PCT / 100.0
        if (!gOn && Float64(T) > Float64(g.T) * (1.0 + x) && T - g.T >= GUARD_ABS) {
            gOn = true
        } else if (gOn && Float64(T) < Float64(g.T) * (1.0 + x / 2.0)) {
            gOn = false
        }
        return gOn
    }
    // чистый best fit: самый короткий подходящий промежуток (среди равных — левый); иначе хвост с минимальной докопкой
    func placePure(L: Int64): Unit {   // результат — zN, zD
        if (F < T) { F = T }   // потоки в этом режиме не работают: граничный промежуток становится обычным
        var ts = tailStart()
        var bgs: Int64 = -1
        var bs: Int64 = 0
        var s = L
        while (s < NB && bgs < 0) {
            let g = firstGap(s, ts)
            if (g >= 0) { bgs = g; bs = s }
            s += 1
        }
        if (bgs < 0) {
            let g = firstBig(ts)
            if (g >= 0) { bgs = g; bs = outLen }
        }
        zN = 0
        if (bgs >= 0) {
            if (bs > L) {
                let lp = pidLeftOf(bgs)
                let rp = pidRightOf(bgs + bs)
                if (lp > 0 && rp > 0) { raceAdd(lp, rp) }
            }
            zD = chooseSide(bgs, bs, L, -1)
            takeFrom(bgs, zD, L)
        } else {
            let tl = T - ts
            if (tl < L) { zN = (L - tl + 19) / 20; dig(zN); ts = tailStart() }
            zD = ts
            takeFrom(ts, zD, L)
        }
    }

    // ---- команды; результат: out[0] = N, out[1] = D
    func place(L: Int64, out: Array<Int64>): Unit {
        shadowCheck()
        if (let Some(sh) <- shadow) { sh.place(L, out) }
        if (!pureBF && !gDropped && gsh.isNone() && npid == 1) {
            let g = Alloc()
            g.pureBF = true
            g.noStreams = true
            gsh = Some(g)
        }
        if (let Some(g) <- gsh) { g.place(L, out) }
        occ += L
        let m = minRecent()
        let pid = npid
        if (L == lastL) {
            runIdx += 1
        } else {
            if (runIdx + 1 >= FIT_RUN) { runPrev = runLast; runLast = runIdx + 1 }
            runIdx = 0
        }
        lastL = L
        feat[pid] = (runIdx & 1) | ((pid & 1) << 1) | (1 << (2 + runIdx % 3)) | (1 << (5 + runIdx % 4)) |
            (1 << (9 + runIdx % 5))
        let bf = pureBF || guardOn()
        let st = if (bf) { -1 } else { chooseStream(pid) }
        stream[pid] = st
        if (st >= 0) { bornLab[st] += 1.0 }
        for (f in 0..NF) {
            birthV[f * 2] *= 0.999; birthV[f * 2 + 1] *= 0.999
            deathV[f * 2] *= 0.999; deathV[f * 2 + 1] *= 0.999
            birthV[f * 2 + ((feat[pid] >> f) & 1)] += 1.0
        }
        if (F > T) { F = T }
        var ts = tailStart()
        let selfHz = if (recDelTot >= 1.0) { rateOfLen(L) } else { 0.0 }
        var N: Int64 = 0
        var D: Int64 = -1
        if (bf) { placePure(L); finishPlace(pid, L, zN, zD, out); return }
        if (zoneTry(L, st, m)) { finishPlace(pid, L, zN, zD, out); return }

        verB = bornLab[1] >= 20.0 && bornLab[0] >= 20.0 && diedLab[1] / bornLab[1] - diedLab[0] / bornLab[0] > 0.9

        var bestC: Float64 = 1e18
        var bgs: Int64 = -1
        var bs: Int64 = 0
        let multi = st < 0 && T >= MULTI_T
        for (s in L..NB) {
            let g = firstGap(s, ts)
            if (g < 0) { continue }
            let c = if (zoneBlocked(g, s, st)) { 1e17 } else { gapCost(g, s, L, st, m, selfHz) }
            if (c < bestC) { bestC = c; bgs = g; bs = s }
            if (multi && c < 1e17 && c - nbW(s - L) < bestC) {
                let ng = someGaps(s, ts, MULTI_K, gk)
                for (i in 1..ng) {
                    let c2 = if (zoneBlocked(gk[i], s, st)) { 1e17 } else { gapCost(gk[i], s, L, st, m, selfHz) }
                    if (c2 < bestC) { bestC = c2; bgs = gk[i]; bs = s }
                }
            }
        }
        let gb = firstBig(ts)
        if (gb >= 0) {
            let bl = outLen
            let c = if (zoneBlocked(gb, bl, st)) { 1e17 } else { gapCost(gb, bl, L, st, m, selfHz) }
            if (c < bestC) { bestC = c; bgs = gb; bs = bl }
        }

        var useGap = bgs >= 0 && bestC < 1e17
        if (useGap && st < 0 && T - tailStart() >= L && LW[13] < bestC) { useGap = false }

        if (useGap) {
            if (bs > L) {
                let lp = pidLeftOf(bgs)
                let rp = pidRightOf(bgs + bs)
                if (lp > 0 && rp > 0) { raceAdd(lp, rp) }
            }
            D = chooseSide(bgs, bs, L, st)
            takeFrom(bgs, D, L)
        } else {
            let fgs = if (F > 0 && tag.get(F) < 0) { -tag.get(F) - 1 } else { F }
            let fl = F - fgs
            var chunk = CHUNK
            let sameGroups = runPrev > 0 && runIdx >= FIT_RUN &&
                50.0 * dabs(Float64(runPrev - runLast)) <= Float64(imax(runPrev, runLast))
            if (sameGroups) { chunk = CHUNK_SAME }
            if (!anyDel) { chunk = EARLY_CHUNK }
            if (let Some(sh) <- shadow) {
                var c = CHUNK_MIN
                let sv = sh.freeIn() - freeIn()
                if (CHUNK_PER_L * L > c) { c = CHUNK_PER_L * L }
                if (npid - shadowBorn > CHUNK_LOW_AGE && sv < c && CHUNK_LOW + L < c) { c = CHUNK_LOW + L }
                if (sv > c) { c = sv }
                if (c < chunk) { chunk = c }
                if (sameGroups) {
                    var rem = (runLast - runIdx) * L
                    if (rem < 4 * L) { rem = 4 * L }
                    if (rem < chunk) { chunk = rem }
                }
            }
            var far = st >= 0 && ((st == 1) != flip)
            if (fl >= L) {
                D = if (far) { F - L } else { fgs }
                takeFrom(fgs, D, L)
                if (far) { F = D }
            } else if (F == T) {
                var need = L - fl
                if (far && chunk - fl > need) { need = chunk - fl }
                N = (need + 19) / 20
                if (far) { N = fitChunk(fl, N, L) }
                dig(N)
                D = if (far) { T - L } else { fgs }
                takeFrom(fgs, D, L)
                F = if (far) { D } else { T }
            } else {
                if (st >= 0) { flip = !flip; far = (st == 1) != flip }
                ts = tailStart()
                let tl = T - ts
                var need = L - tl
                if (st >= 0 && chunk - tl > need) { need = chunk - tl }
                N = if (need > 0) { (need + 19) / 20 } else { 0 }
                if (st >= 0 && N > 0) { N = fitChunk(tl, N, L) }
                dig(N)
                ts = tailStart()
                D = if (far) { T - L } else { ts }
                takeFrom(ts, D, L)
                F = if (far) { D } else { T }
            }
        }
        finishPlace(pid, L, N, D, out)
    }
    func finishPlace(pid: Int64, L: Int64, N: Int64, D: Int64, out: Array<Int64>): Unit {
        for (i in 0..L) { owner.set(D + i, pid) }
        pos[pid] = D
        len[pid] = L
        alive[pid] = true
        npid += 1
        bitAdd(pid, 1)
        birthL[L] += 1.0
        aliveCntL[L] += 1
        aliveTot += 1
        lastOfLen[L * 8 + lastPtr[L]] = pid
        lastPtr[L] = (lastPtr[L] + 1) & 7
        lastSeen[L] = pid
        for (l in 1..=MAXL) { cnt[l] *= 0.999; cntF[l] *= 0.97; cntR[l] *= 0.999 }
        cntTot = cntTot * 0.999 + 1.0
        cnt[L] += 1.0
        cntFTot = cntFTot * 0.97 + 1.0
        cntF[L] += 1.0
        cntR[L] += 1.0
        if (cntFTot >= 25.0 && pid - lastShift > SHIFT_GAP) {
            var tv: Float64 = 0.0
            for (l in 1..=MAXL) { tv += dabs(cntF[l] / cntFTot - cntR[l] / cntTot) }
            if (tv / 2.0 > SHIFT_TV) {
                for (l in 1..=MAXL) { cntR[l] = cntF[l] / cntFTot * cntTot }
                lastShift = pid
            }
        }
        lastWasDel = false
        out[0] = N
        out[1] = D
    }

    func remove(pid: Int64): Int64 {
        anyDel = true
        let d = pos[pid]
        let l = len[pid]
        occ -= l
        if (let Some(sh) <- shadow) { sh.remove(pid) }
        if (let Some(g) <- gsh) { g.remove(pid) }
        deathL[l] += 1.0
        if (aliveTot > 1) {
            let q = Float64(aliveTot - bitSum(pid)) / Float64(aliveTot - 1)
            qSum = qSum * 0.99 + q
            qTot = qTot * 0.99 + 1.0
        }
        bitAdd(pid, -1)
        alive[pid] = false
        raceOnDeath(pid)
        let b = ageBucket(npid - 1 - pid)
        deathAge[b] += 1.0
        for (k in 0..=b) { deathSuf[k] += 1.0 }
        for (k in 1..=MAXL) { recDel[k] *= 0.98 }
        recDelTot = recDelTot * 0.98 + 1.0
        recDel[l] += 1.0
        aliveCntL[l] -= 1
        aliveTot -= 1
        var hit = false
        for (k in 0..8) { if (lastOfLen[l * 8 + k] == pid) { hit = true } }
        lifoHit = lifoHit * 0.99 + (if (hit) { 1.0 } else { 0.0 })
        lifoTot = lifoTot * 0.99 + 1.0
        if (zl.get(d) == l) { slots[l].push(d) }
        let e = freeRange(d, l)
        if (d == F) { F = e }
        if (stream[pid] >= 0) { diedLab[stream[pid]] += 1.0 }
        for (f in 0..NF) {
            let v = (feat[pid] >> f) & 1
            deathV[f * 2 + v] += 1.0
            if (lastWasDel) {
                coSame[f] *= 0.998
                coTot[f] *= 0.998
                coTot[f] += 1.0
                coSame[f] += (if (v == ((feat[lastDelPid] >> f) & 1)) { 1.0 } else { 0.0 })
            }
        }
        lastWasDel = true
        lastDelPid = pid
        return d
    }
}

// ---------------------------------------------------------------- ввод-вывод
// Без String на команду (readln/Int64.parse/интерполяция строк дают мусор на каждую из 200 000 команд):
// число читается по Rune из ConsoleReader (внутри он читает ровно строку и не ждёт больше — интерактив не блокируется),
// ответ собирается в байтовом буфере и уходит одним write + flush.
let obuf: Array<Byte> = Array<Byte>(64, repeat: 0)
let stdin = getStdIn()
let stdout = getStdOut()

// -2^63 — конец ввода
func readInt(): Int64 {
    var c: Int64 = 0
    while (true) {
        match (stdin.read()) {
            case Some(r) => c = Int64(UInt32(r))
            case None => return Int64.Min
        }
        if (c == 45 || (c >= 48 && c <= 57)) { break }
    }
    var neg = false
    if (c == 45) {
        neg = true
        match (stdin.read()) {
            case Some(r) => c = Int64(UInt32(r))
            case None => return Int64.Min
        }
    }
    var x: Int64 = 0
    while (c >= 48 && c <= 57) {
        x = x * 10 + (c - 48)
        match (stdin.read()) {
            case Some(r) => c = Int64(UInt32(r))
            case None => c = 10
        }
    }
    return if (neg) { -x } else { x }
}

func putNum(p: Int64, v0: Int64): Int64 {
    var v = v0
    var n = p
    if (v == 0) {
        obuf[n] = 48
        return n + 1
    }
    while (v > 0) {
        obuf[n] = UInt8(48 + v % 10)
        n += 1
        v /= 10
    }
    var i = p
    var j = n - 1
    while (i < j) {
        let t = obuf[i]
        obuf[i] = obuf[j]
        obuf[j] = t
        i += 1
        j -= 1
    }
    return n
}

main(): Int64 {
    let A = Alloc()
    let res = Array<Int64>(2, repeat: 0)
    while (true) {
        let c = readInt()
        if (c == Int64.Min || c == 0) { return 0 }   // конец: ничего не печатаем
        var n: Int64 = 0
        if (c > 0) {
            A.place(c, res)
            n = putNum(0, res[0])
            obuf[n] = 32
            n = putNum(n + 1, res[1])
        } else {
            n = putNum(0, A.remove(-c))
        }
        obuf[n] = 10
        n += 1
        stdout.write(obuf[0..n])
        stdout.flush()             // обязателен после каждого ответа
    }
    return 0
}
