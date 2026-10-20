// scorefit — порт solutions/scorefit.cpp на Cangjie (仓颉) 1.0.
//
// Проверено на cjc 1.0.5 (2026-10-20): компилируется без предупреждений, ответы побитово = C++-версии из
// коммита 5d99f4b (порт соответствует ей, а НЕ текущей solutions/scorefit.cpp) на 46 тестах, через интерактор
// 0 WA/TL; худший случай (200 тыс. труб по 40 м) — 0.5 с, 113 МБ.
// Сборка на macOS: `source build/cangjie-sdk/cangjie/envsetup.sh`, затем
//   SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX15.4.sdk cjc -O2 scorefit.cj -o scorefit
// (lld из SDK не читает MacOSX26.sdk). Сравнить ответы с C++-версией:
//   for f in materials/tests/*.in; do tail -n +2 $f > /tmp/c.txt; ./scorefit < /tmp/c.txt > /tmp/a.txt;
//     build/bin/scorefit < /tmp/c.txt > /tmp/b.txt; cmp /tmp/a.txt /tmp/b.txt || echo DIFF $f; done
// Ответы должны совпадать побитово (логика и порядок операций с плавающей точкой — как в C++).
//
// Отличия от C++ только технические: все целые — Int64 (страницы позиций — Int32); таблица страниц
// фиксированного размера (хватает на траншею до 16.7 млн м, максимум по задаче ≈ 8.5 млн); массивы по PID —
// фиксированного размера 200 002 (команд ≤ 200 000).
package scorefit

import std.env.*
import std.convert.*

let NB: Int64 = 128          // промежутки короче NB — в корзинах по точной длине
let MAXL: Int64 = 40
let CHUNK: Int64 = 400
let MAXPID: Int64 = 200002
let PSH: Int64 = 16          // размер страницы 2^16
let PMASK: Int64 = 65535
let NPAGES: Int64 = 256

func dmax(a: Float64, b: Float64): Float64 { if (a > b) { a } else { b } }
func dmin(a: Float64, b: Float64): Float64 { if (a < b) { a } else { b } }
func dabs(a: Float64): Float64 { if (a < 0.0) { -a } else { a } }
func imax(a: Int64, b: Int64): Int64 { if (a > b) { a } else { b } }

// ---------------------------------------------------------------- постраничный массив (Int32 — ради памяти)
class Paged {
    var pg: Array<Array<Int32>> = Array<Array<Int32>>(NPAGES, repeat: Array<Int32>())
    var npg: Int64 = 0
    func grow(n: Int64): Unit {
        while ((npg << PSH) < n) {
            pg[npg] = Array<Int32>(PMASK + 1, repeat: 0)
            npg += 1
        }
    }
    func get(i: Int64): Int64 { Int64(pg[i >> PSH][i & PMASK]) }
    func set(i: Int64, v: Int64): Unit { pg[i >> PSH][i & PMASK] = Int32(v) }
}

// ---------------------------------------------------------------- двоичная куча (минимум)
class Heap {
    var a: Array<Int64> = Array<Int64>(16, repeat: 0)
    var n: Int64 = 0
    func empty(): Bool { n == 0 }
    func top(): Int64 { a[0] }
    func push(v: Int64): Unit {
        if (n == a.size) {
            let b = Array<Int64>(a.size * 2, repeat: 0)
            for (k in 0..n) { b[k] = a[k] }
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
    let tag = Paged()     // >0: начало промежутка длины tag; <0: конец промежутка, начавшегося в -tag-1
    let owner = Paged()   // PID трубы в клетке
    let bucket: Array<Heap> = Array<Heap>(NB, {_ => Heap()})
    let big = Heap()

    let pos = Array<Int64>(MAXPID, repeat: 0)
    let len = Array<Int64>(MAXPID, repeat: 0)
    let feat = Array<Int64>(MAXPID, repeat: 0)
    let stream = Array<Int64>(MAXPID, repeat: -1)
    var npid: Int64 = 1   // следующий PID

    let cnt = Array<Float64>(MAXL + 1, repeat: 0.0)
    var cntTot: Float64 = 0.0
    let birthL = Array<Float64>(MAXL + 1, repeat: 0.0)
    let deathL = Array<Float64>(MAXL + 1, repeat: 0.0)
    let recDel = Array<Float64>(MAXL + 1, repeat: 0.0)
    var recDelTot: Float64 = 0.0
    let aliveCntL = Array<Int64>(MAXL + 1, repeat: 0)
    var aliveTot: Int64 = 0

    var F: Int64 = 0
    var runIdx: Int64 = 0
    var lastL: Int64 = 0
    let NF: Int64 = 5
    let coSame = Array<Float64>(5, repeat: 0.0)
    let coTot = Array<Float64>(5, repeat: 0.0)
    let birthV = Array<Float64>(10, repeat: 0.0)   // [f*2 + v]
    let deathV = Array<Float64>(10, repeat: 0.0)
    var curFeat: Int64 = -1
    let bornLab = Array<Float64>(2, repeat: 0.0)
    let diedLab = Array<Float64>(2, repeat: 0.0)
    var lastDelPid: Int64 = 0
    var lastWasDel: Bool = false

    init() { ensure(1024) }

    // ---- промежутки
    func ensure(n: Int64): Unit { tag.grow(n + 1); owner.grow(n + 1) }
    func addGap(s: Int64, l: Int64): Unit {
        tag.set(s, l); tag.set(s + l, -s - 1)
        if (l < NB) { bucket[l].push(s) } else { big.push((l << 32) | s) }
    }
    func delGap(s: Int64): Unit {
        let l = tag.get(s)
        tag.set(s, 0); tag.set(s + l, 0)
    }
    func freeRange(s0: Int64, l: Int64): Int64 {
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
        let h = bucket[s]
        let keep = Array<Int64>(3, repeat: 0)
        var nk: Int64 = 0
        var res: Int64 = -1
        while (!h.empty() && nk < 3) {
            let g = h.top()
            if (tag.get(g) != s || (nk > 0 && keep[nk - 1] == g)) { h.pop(); continue }
            h.pop()
            keep[nk] = g
            nk += 1
            if (g != ts && g + s != F) { res = g; break }
        }
        for (i in 0..nk) { h.push(keep[i]) }
        return res
    }
    // возвращает start или -1; длина — во втором элементе out
    func firstBig(ts: Int64, out: Array<Int64>): Int64 {
        let keep = Array<Int64>(3, repeat: 0)
        var nk: Int64 = 0
        var res: Int64 = -1
        while (!big.empty() && nk < 3) {
            let k = big.top()
            let l = k >> 32
            let g = k & 0xFFFFFFFF
            if (tag.get(g) != l || (nk > 0 && keep[nk - 1] == k)) { big.pop(); continue }
            big.pop()
            keep[nk] = k
            nk += 1
            if (g != ts && g + l != F) { res = g; out[0] = l; break }
        }
        for (i in 0..nk) { big.push(keep[i]) }
        return res
    }

    // ---- статистика
    func minRecent(): Int64 {
        for (l in 1..=MAXL) { if (cnt[l] > 0.002 * cntTot) { return l } }
        return 1
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

    // ---- потоки
    func kappa(f: Int64): Float64 {
        if (coTot[f] < 20.0) { return 0.0 }
        let b0 = birthV[f * 2]
        let b1 = birthV[f * 2 + 1]
        let bt = b0 + b1
        if (bt <= 0.0) { return 0.0 }
        let mn = if (b0 < b1) { b0 } else { b1 }
        if (mn < 0.3 * bt) { return 0.0 }
        let e = (b0 * b0 + b1 * b1) / (bt * bt)
        let o = coSame[f] / coTot[f]
        return if (e >= 1.0) { 0.0 } else { (o - e) / (1.0 - e) }
    }
    func pdie(f: Int64, v: Int64): Float64 { (deathV[f * 2 + v] + 1.0) / (birthV[f * 2 + v] + 2.0) }
    func chooseStream(pid: Int64): Int64 {
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
        return if (lp > rp) { gs + r } else { gs }
    }

    func gapCost(gs: Int64, s: Int64, L: Int64, st: Int64, m: Int64, selfHz: Float64): Float64 {
        if (st < 0 && s > L && s - L < m && s - L >= 4) { return 1e17 }
        if (T < 60) {
            if (s == L) { return -1000.0 + 2000.0 * dmax(pdieLen(pidLeftOf(gs)), pdieLen(pidRightOf(gs + s))) }
            return -Float64(s)
        }
        if (s != L) { return Float64(s - L) }
        var h: Float64 = 0.0
        let nb0 = pidLeftOf(gs)
        let nb1 = pidRightOf(gs + s)
        for (k in 0..2) {
            let nb = if (k == 0) { nb0 } else { nb1 }
            if (nb <= 0) { continue }
            let mg = s + len[nb]
            if (mg > MAXL || cnt[mg] < 0.05 * cntTot) { continue }
            h = dmax(h, hzRate(nb))
        }
        h -= selfHz
        return if (h >= 0.5) { 3.0 * h } else { 0.0 }
    }

    // ---- команды; результат: out[0] = N, out[1] = D
    func place(L: Int64, out: Array<Int64>): Unit {
        let m = minRecent()
        let pid = npid
        npid += 1
        if (L == lastL) { runIdx += 1 } else { runIdx = 0 }
        lastL = L
        feat[pid] = (runIdx & 1) | ((pid & 1) << 1) | (1 << (2 + runIdx % 3))
        let st = chooseStream(pid)
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

        var bestC: Float64 = 1e18
        var bgs: Int64 = -1
        var bs: Int64 = 0
        for (s in L..NB) {
            let g = firstGap(s, ts)
            if (g < 0) { continue }
            let c = gapCost(g, s, L, st, m, selfHz)
            if (c < bestC) { bestC = c; bgs = g; bs = s }
        }
        let bo = Array<Int64>(1, repeat: 0)
        let gb = firstBig(ts, bo)
        if (gb >= 0) {
            let c = gapCost(gb, bo[0], L, st, m, selfHz)
            if (c < bestC) { bestC = c; bgs = gb; bs = bo[0] }
        }

        var useGap = bgs >= 0 && bestC < 1e17
        if (st == 1 && bornLab[1] >= 20.0 && bornLab[0] >= 20.0 &&
            diedLab[1] / bornLab[1] - diedLab[0] / bornLab[0] > 0.5) { useGap = false }

        var N: Int64 = 0
        var D: Int64 = 0
        if (useGap) {
            D = chooseSide(bgs, bs, L, st)
            takeFrom(bgs, D, L)
        } else {
            let fgs = if (F > 0 && tag.get(F) < 0) { -tag.get(F) - 1 } else { F }
            let fl = F - fgs
            if (fl >= L) {
                D = if (st == 1) { F - L } else { fgs }
                takeFrom(fgs, D, L)
                if (st == 1) { F = D }
            } else if (F == T) {
                var need = L - fl
                if (st == 1 && CHUNK - fl > need) { need = CHUNK - fl }
                N = (need + 19) / 20
                dig(N)
                D = if (st == 1) { T - L } else { fgs }
                takeFrom(fgs, D, L)
                F = if (st == 1) { D } else { T }
            } else {
                ts = tailStart()
                let tl = T - ts
                var need = L - tl
                if (st >= 0 && CHUNK - tl > need) { need = CHUNK - tl }
                N = if (need > 0) { (need + 19) / 20 } else { 0 }
                dig(N)
                ts = tailStart()
                D = if (st == 1) { T - L } else { ts }
                takeFrom(ts, D, L)
                F = if (st == 1) { D } else { T }
            }
        }
        for (i in 0..L) { owner.set(D + i, pid) }
        pos[pid] = D
        len[pid] = L
        birthL[L] += 1.0
        aliveCntL[L] += 1
        aliveTot += 1
        for (l in 1..=MAXL) { cnt[l] *= 0.999 }
        cntTot = cntTot * 0.999 + 1.0
        cnt[L] += 1.0
        lastWasDel = false
        out[0] = N
        out[1] = D
    }

    func remove(pid: Int64): Int64 {
        let d = pos[pid]
        let l = len[pid]
        deathL[l] += 1.0
        for (k in 1..=MAXL) { recDel[k] *= 0.98 }
        recDelTot = recDelTot * 0.98 + 1.0
        recDel[l] += 1.0
        aliveCntL[l] -= 1
        aliveTot -= 1
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
                if (v == ((feat[lastDelPid] >> f) & 1)) { coSame[f] += 1.0 }
            }
        }
        lastWasDel = true
        lastDelPid = pid
        return d
    }
}

main(): Int64 {
    let input = getStdIn()
    let output = getStdOut()
    let A = Alloc()
    let res = Array<Int64>(2, repeat: 0)
    while (true) {
        // каждая команда интерактора — отдельная строка
        let line = input.readln()
        match (line) {
            case Some(text) =>
                if (text.size == 0) { continue }
                let c = Int64.parse(text)
                if (c == 0) { return 0 }   // конец: ничего не печатаем
                if (c > 0) {
                    A.place(c, res)
                    output.write("${res[0]} ${res[1]}\n")
                } else {
                    output.write("${A.remove(-c)}\n")
                }
                output.flush()             // обязателен после каждого ответа
            case None => return 0
        }
    }
    return 0
}
