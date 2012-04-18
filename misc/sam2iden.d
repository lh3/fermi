#!/usr/bin/env rdmd

import std.stdio, std.array, bio, std.conv;

void main(string[] args)
{
	if (args.length == 1) {
		writeln("\nUsage: sam2iden.d <in.sam>\n");
		writeln("Output format: queryName, queryStart, queryEnd, strand, refStart, refEnd, mapQ, blastIdentity, blatIdentity\n");
		return;
	}
	auto fp = new ZFile(args[1]);
	ubyte[] data;
	while (fp.readto(data) > 0) {
		if (data[0] == '@') continue;
		auto t = split(cast(string)data);
		auto flag = to!int(t[1]);
		if (flag & 4) continue;
		auto cs = sam_parse_cigar(t[5]);
		int ndiff;
		for (int i = 11; i < cast(int)t.length; ++i)
			if (t[i][0..5] == "NM:i:")
				ndiff = to!int(t[i][5..$]);
		auto qlen = cs.n_M + cs.clip[0] + cs.clip[1] + cs.n_I; // qlen is the query sequence length
		auto pos = to!int(t[3]) - 1;
		write(t[0], '\t');
		if (flag&16) write(cs.clip[1], '\t', qlen - cs.clip[0], "\t-\t");
		else write(cs.clip[0], '\t', qlen - cs.clip[1], "\t+\t");
		qlen -= cs.clip[0] + cs.clip[1]; // qlen is now the alignment length
		write(t[2], '\t', pos, '\t', pos + cs.n_M + cs.n_D, '\t', t[4], '\t',
			  cast(double)(qlen + cs.n_D - ndiff) / (qlen + cs.n_D), '\t',
			  cast(double)(qlen - cs.n_I - (ndiff - cs.n_I - cs.n_D)) / (qlen - cs.n_I), '\n');
	}
}
