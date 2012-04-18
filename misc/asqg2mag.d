#!/usr/bin/env rdmd

import std.array, std.stdio, std.conv, bio;

void asqg2mag(string fn) {
	struct ipair_t {
		ulong idd;
		int o;
		this(ulong x, int y) { idd = x; o = y; }
	}
	struct asqgv_t {
		string s;
		ipair_t[][2] nei;
		this(string t) { s = t; }
	}

	ulong[string] v;
	asqgv_t[] seqs;
	ubyte[] line;
	auto f = new ZFile(fn);
	while (f.readto(line) >= 0) {
		auto t = split(cast(string)line);
		if (t[0] == "VT") {
			v[t[1]] = seqs.length;
			seqs ~= asqgv_t(t[2]);
		} else if (t[0] == "ED") {
			int[6] x; // start1, end1, len1, start2, end2, len2
			for (int i = 3; i <= 8; ++i)
				x[i - 3] = to!int(t[i]);
			++x[1]; ++x[4];
			int o = x[1] - x[0]; // overlap length
			assert(o == x[4] - x[3]); // no gaps are allowed
			ulong id1 = v[t[1]], id2 = v[t[2]];
			int y1 = x[0] == 0? 0 : x[2] - x[1] == 0? 1 : -1; // which end is linked
			int y2 = x[3] == 0? 0 : x[5] - x[4] == 0? 1 : -1;
			assert(y1 != -1 && y2 != -1); // only works for end-to-end overlap
			seqs[id1].nei[y1] ~= ipair_t(id2<<1|y2, o);
			seqs[id2].nei[y2] ~= ipair_t(id1<<1|y1, o);
		}
	}
	for (ulong i = 0; i < seqs.length; ++i) {
		write(">", i<<1, ":", i<<1|1, "\t1");
		for (int j = 0; j < 2; ++j) {
			write('\t');
			auto p = seqs[i].nei[j];
			if (p.length) {
				for (int k = 0; k < p.length; ++k)
					write(p[k].idd, ',', p[k].o, ';');
			} else write('.');
		}
		write('\n');
		writeln(seqs[i].s);
	}
}

void main(string[] args) {
	if (args.length == 1) {
		writeln("Usage: asqg2mag <graph.asqg.gz>");
		return;
	}
	asqg2mag(args[1]);
}
