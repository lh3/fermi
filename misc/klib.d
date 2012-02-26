module klib;

import std.string, std.ascii, std.c.stdlib, std.c.string;

extern(C) void *gzopen(const char *fn, const char *mode);
extern(C) void *gzdopen(int fd, const char *mode);
extern(C) int gzread(void *fp, void *buf, uint len);
extern(C) int gzclose(void *fp);

class ZFile {
	int _begin, _end;
	bool _eof;
	uint _buf_size;
	void *_fp;
	ubyte *_buf;
public:
	this(string fn = "", int buf_size = 16384) {
		_begin = _end = 0, _eof = false;
		_buf_size = buf_size;
		_buf = cast(ubyte*)std.c.stdlib.malloc(_buf_size);
		_fp = fn.length? gzopen(std.string.toStringz(fn), "rb") : gzdopen(0, "rb");
	}
	~this() {
		std.c.stdlib.free(_buf);
		gzclose(_fp);
	}
	int readto(ref ubyte[] dat, ubyte delimiter = '\n', bool append = false) {
		if (!append) dat.length = 0;
		if (_begin >= _end && _eof) return -1;
		while (1) {
			if (_begin >= _end) {
				if (!_eof) {
					_begin = 0;
					_end = gzread(_fp, cast(void*)_buf, _buf_size);
					if (_end < _buf_size) _eof = true;
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
