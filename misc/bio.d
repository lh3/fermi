module bio;

import std.string, std.ascii, std.c.stdlib, std.c.string;

/*******************************
 * Buffered gzip stream reader *
 *******************************/

extern(C) void *gzopen(const char *fn, const char *mode);
extern(C) void *gzdopen(int fd, const char *mode);
extern(C) int gzread(void *fp, void *buf, uint len);
extern(C) int gzclose(void *fp);

class ZFile {
	int _begin, _end;
	uint _buf_size;
	void *_fp;
	ubyte *_buf;
public:
	bool eof;
	void open(string fn = "-", int buf_size = 16384) {
		if (_fp != null) this.close();
		_begin = _end = 0, eof = false;
		_buf_size = buf_size;
		_buf = cast(ubyte*)std.c.stdlib.malloc(_buf_size);
		_fp = fn.length && fn != "-"? gzopen(std.string.toStringz(fn), "rb") : gzdopen(0, "rb");
		if (_fp == null) eof = true;
	}
	void close() {
		if (_fp != null) {
			std.c.stdlib.free(_buf);
			gzclose(_fp);
			_fp = null;
		}
	}
	this(string fn, int buf_size = 16384) { this.open(fn, buf_size); }
	this() { eof = true; _fp = null; }
	~this() { this.close(); }
	int readto(ref ubyte[] dat, ubyte delimiter = '\n', bool append = false) {
		if (!append) dat.length = 0;
		if (_begin >= _end && eof) return -1;
		while (1) {
			if (_begin >= _end) {
				if (!eof) {
					_begin = 0;
					_end = gzread(_fp, cast(void*)_buf, _buf_size);
					if (_end < _buf_size) eof = true;
					if (_end == 0) break;
				} else break;
			}
			int i;
			if (delimiter > 2) {
				for (i = _begin; i < _end; ++i)
					if (_buf[i] == delimiter) break;
			} else if (delimiter == 0) {
				for (i = _begin; i < _end; ++i)
					if (std.ascii.isWhite(_buf[i])) break;
			} else if (delimiter == 1) {
				for (i = _begin; i < _end; ++i)
					if (std.ascii.isWhite(_buf[i]) && _buf[i] != ' ') break;
			}
			auto old_l = dat.length;
			dat.length += i - _begin;
			std.c.string.memcpy(dat.ptr + old_l, _buf + _begin, i - _begin);
			_begin = i + 1;
			if (i < _end) return cast(int)_buf[i];
		}
		return -1;
	}
}

/*******************************
 * Bio sequence transformation *
 *******************************/

char[128] seq_comp_table = [
      0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
     16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
     32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
     48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
     64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
    'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
     64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
    'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
];

char[] seq_revcomp(char[] seq) {
	char[] rev;
	rev.length = seq.length;
	for (int i = cast(int)seq.length - 1; i >= 0; --i)
		rev[seq.length - 1 - i] = seq_comp_table[cast(int)seq[i]];
	return rev;
}

/*************************
 * SAM related ulitities *
 *************************/

struct sam_cigarsum_t {
	int n_M, n_I, n_D, n_N;
	int[2] clip;
}

sam_cigarsum_t *sam_parse_cigar(string cigar)
{
	auto p = new sam_cigarsum_t;
	auto q = cigar.ptr;
	while (q - cigar.ptr < cigar.length) {
		char *r;
		auto oplen = strtol(q, &r, 10);
		if (*r == 'S' || *r == 'H') {
			p.clip[q == cigar.ptr? 0 : 1] = cast(int)oplen;
		} else if (*r == 'M') p.n_M += oplen;
		else if (*r == 'I') p.n_I += oplen;
		else if (*r == 'D') p.n_D += oplen;
		else if (*r == 'N') p.n_N += oplen;
		q = cast(immutable(char)*)r + 1;
	}
	return p;
}
