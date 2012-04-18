#!/usr/bin/env rdmd

import std.stdio, std.array, std.conv, std.getopt, std.algorithm, std.c.stdio, bio;

struct cmdopt_t {
	bool is_print;
	int min_len, max_gap, min_q;
	double mask_level;
}

struct stats_t {
	ulong L, n_un, l_un, n_dropped;
	int[5] n_b, n_bg;
	int[] len;
}

struct aln_t {
	string sam, chr;
	int pos, len, qlen, rlen, flag, mapq, qbeg;
	int[2] clip;
}

aln_t *parse_aln(string l, ref string[] t)
{
	auto p = new aln_t;
	p.sam = l;
	p.chr = t[2];
	p.pos = to!int(t[3]) - 1;
	p.mapq = to!int(t[4]);
	p.flag = to!int(t[1]);
	if ((p.flag&4) == 0) { // mapped
		auto cs = sam_parse_cigar(t[5]);
		p.qlen = cs.n_M + cs.n_I;
		p.rlen = cs.n_M + cs.n_D + cs.n_N;
		p.clip[0] = cs.clip[0]; p.clip[1] = cs.clip[1];
		p.qbeg = p.clip[!!(p.flag&16)];
		p.len = p.clip[0] + p.clip[1] + p.qlen;
	} else {
		auto s = split(l);
		p.len = cast(int)s[9].length;
	}
	return p;
}

void count_break(ref int[5] c, ref aln_t*[] a, const ref cmdopt_t opt)
{
	int[5] b = [cast(int)a.length, 0, 0, 0, 0];
	foreach (p; a) {
		if (p.mapq < opt.min_q) continue;
		++b[1];
		if (p.qlen >= 100) {
			++b[2];
			if (p.qlen >= 200) {
				++b[3];
				if (p.qlen >= 500) ++b[4];
			}
		}
	}
	for (int i = 0; i < 5; ++i)
		if (b[i]) c[i] += b[i] - 1;
}

void analyze_aln(ref aln_t*[] a, ref stats_t s, const ref cmdopt_t opt)
{
	// special treatment of unmapped
	if (a.length == 1 && (a[0].flag&4) == 4) {
		++s.n_un; s.l_un += a[0].len;
		if (opt.is_print) writeln(a[0].sam);
		return;
	}
	// apply mask_level
	if (a.length > 1) { // multi-part alignment
		aln_t*[] tmp;
		foreach (p; a) {
			bool dropped = false;
			foreach (q; tmp) {
				auto beg = p.qbeg > q.qbeg? p.qbeg : q.qbeg;
				auto end = p.qbeg+p.qlen < q.qbeg+q.qlen? p.qbeg+p.qlen : q.qbeg+q.qlen;
				if (beg < end && cast(double)(end - beg) > p.qlen * opt.mask_level) {
					dropped = true;
					break;
				}
			}
			if (!dropped) tmp ~= p;
			else ++s.n_dropped;
		}
		a = tmp;
		count_break(s.n_b, a, opt);
	}
	foreach (p; a) s.len ~= p.qlen;
	if (opt.is_print)
		foreach (p; a) writeln(p.sam);
	// patch small gaps
	if (a.length > 1) { // still multi-part
		sort!((x,y){return x.chr < y.chr || (x.chr == y.chr && x.pos < y.pos);})(a);
		for (int i = 1; i < a.length; ++i) {
			auto p = a[i], q = a[i-1];
			if (p.chr == q.chr && (p.flag&16) == (q.flag&16)) { // same chr and same flag
				auto gapr = p.pos - (q.pos + q.rlen);
				auto gapq = p.clip[0] - (q.clip[0] + q.qlen);
				if (gapr < 0) gapr = -gapr;
				if (gapq < 0) gapq = -gapq;
				if (gapr < opt.max_gap && gapq < opt.max_gap) {
					p.qlen = p.clip[0] + p.qlen - q.clip[0]; p.clip[0] = q.clip[0];
					p.rlen = p.pos + p.rlen - q.pos; p.pos = q.pos;
					q.flag |= 4; // delete a[i-1]
				}
			}
		}
		aln_t*[] tmp;
		foreach (p; a)
			if ((p.flag&4) == 0) tmp ~= p;
		a = tmp;
		count_break(s.n_bg, a, opt);
	}
}

void main(string[] args)
{
	cmdopt_t opt = {false, 150, 500, 10, 0.5};
	getopt(args, std.getopt.config.bundling, "l", &opt.min_len, "q", &opt.min_q, "p", &opt.is_print, "m", &opt.mask_level);
	if (args.length == 1) {
		writeln("\nUsage:   sam2break.d <contig-aln.sam>\n");
		writeln("Options: -l INT     exclude contigs shorter than INT bp [", opt.min_len, ']');
		writeln("         -q INT     exclude alignments with maqQ lower than INT [", opt.min_q, ']');
		writeln("         -m FLOAT   exclude alignments overlapping with a long alignment by FLOAT fraction [", opt.mask_level, ']');
		writeln("         -g INT     join alignments separated by a gap shorter than INT bp [", opt.max_gap, "]\n");
		return;
	}
	auto f = new ZFile(args[1]);
	ubyte[] l;
	string last;
	aln_t*[] a;
	stats_t s;
	while (f.readto(l) >= 0) {
		if (l[0] == '@') {
			if (opt.is_print) writeln(cast(string)l);
			continue;
		}
		int col = 0, l_end;
		for (l_end = 0; l_end < l.length; ++l_end)
			if (l[l_end] == '\t' && (++col) == 6) break;
		auto t = split(cast(string)l[0..l_end]); // do not split beyond the 6th TAB; this is for efficiency
		if (t[0] != last) {
			analyze_aln(a, s, opt);
			a.length = 0;
			last = t[0];
		}
		auto p = parse_aln(cast(string)l, t);
		if (p.len >= opt.min_len) a ~= p;
	}
	analyze_aln(a, s, opt);
	if (!opt.is_print) {
		ulong L, N50, len;
		sort!((a, b){return a > b;})(s.len);
		foreach (p; s.len) L += p;
		foreach (p; s.len)
			if ((len += p) >= L/2) {
				N50 = p; break;
			}
		writeln("Number of unmapped contigs: ", s.n_un);
		writeln("Total length of unmapped contigs: ", s.l_un);
		writeln("Number of alignments dropped due to excessive overlaps: ", s.n_dropped);
		writeln("Mapped contig bases: ", L);
		writeln("Mapped N50: ", N50);
		writeln("Number of break points: ", s.n_b[0]);
		writefln("Number of Q%d break points longer than (0,100,200,500)bp: (%d,%d,%d,%d)", opt.min_q, s.n_b[1], s.n_b[2], s.n_b[3], s.n_b[4]);
		writefln("Number of break points after patching gaps short than %dbp: %d", opt.max_gap, s.n_bg[0]);
		writefln("Number of Q%d break points longer than (0,100,200,500)bp after gap patching: (%d,%d,%d,%d)", opt.min_q, s.n_bg[1], s.n_bg[2], s.n_bg[3], s.n_bg[4]);
	}
}
