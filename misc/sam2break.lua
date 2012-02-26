#!/usr/bin/env luajit

--[[ klib.lua routines ]]--

function os.getopt(args, ostr)
	local arg, place = nil, 0;
	return function ()
		if place == 0 then -- update scanning pointer
			place = 1
			if #args == 0 or args[1]:sub(1, 1) ~= '-' then place = 0; return nil end
			if #args[1] >= 2 then
				place = place + 1
				if args[1]:sub(2, 2) == '-' then -- found "--"
					place = 0
					table.remove(args, 1);
					return nil;
				end
			end
		end
		local optopt = args[1]:sub(place, place);
		place = place + 1;
		local oli = ostr:find(optopt);
		if optopt == ':' or oli == nil then -- unknown option
			if optopt == '-' then return nil end
			if place > #args[1] then
				table.remove(args, 1);
				place = 0;
			end
			return '?';
		end
		oli = oli + 1;
		if ostr:sub(oli, oli) ~= ':' then -- do not need argument
			arg = nil;
			if place > #args[1] then
				table.remove(args, 1);
				place = 0;
			end
		else -- need an argument
			if place <= #args[1] then  -- no white space
				arg = args[1]:sub(place);
			else
				table.remove(args, 1);
				if #args == 0 then -- an option requiring argument is the last one
					place = 0;
					if ostr:sub(1, 1) == ':' then return ':' end
					return '?';
				else arg = args[1] end
			end
			table.remove(args, 1);
			place = 0;
		end
		return optopt, arg;
	end
end

function string:split(sep, n)
	local a, start = {}, 1;
	sep = sep or "%s+";
	repeat
		local b, e = self:find(sep, start);
		if b == nil then
			table.insert(a, self:sub(start));
			break
		end
		a[#a+1] = self:sub(start, b - 1);
		start = e + 1;
		if n and #a == n then
			table.insert(a, self:sub(start));
			break
		end
	until start > #self;
	return a;
end

function io.xopen(fn, mode)
	mode = mode or 'r';
	if fn == nil then return io.stdin;
	elseif fn == '-' then return (mode == 'r' and io.stdin) or io.stdout;
	elseif fn:sub(-3) == '.gz' then return (mode == 'r' and io.popen('gzip -dc ' .. fn, 'r')) or io.popen('gzip > ' .. fn, 'w');
	elseif fn:sub(-4) == '.bz2' then return (mode == 'r' and io.popen('bzip2 -dc ' .. fn, 'r')) or io.popen('bgzip2 > ' .. fn, 'w');
	else return io.open(fn, mode) end
end

--[[ ]]--

local bit = require('bit')

function analyze(stack, stats, opt)

	function count_break(stack, opt)
		local b = {#stack, 0, 0, 0, 0}
		for _, v in ipairs(stack) do
			if v[5] >= opt.min_q then
				local l = v[9] - v[8]
				b[2] = b[2] + 1
				if l >= 100 then
					b[3] = b[3] + 1
					if l >= 200 then
						b[4] = b[4] + 1
						if l >= 500 then
							b[5] = b[5] + 1
						end
					end
				end
			end
		end
		return b
	end

	-- special treatment of unmapped reads
	if #stack == 1 and bit.band(stack[1][2], 4) == 4 then
		stats.n_un, stats.l_un = stats.n_un + 1, stats.l_un + stack[1][7]
		if opt.is_print then print(stack[1].l) end
		return
	end
	-- compute the start and end of the alignment
	for i = 1, #stack do
		local qb, ql, tl = -1, 0, 0
		local cigar = stack[i][6]
		for l in cigar:gmatch('(%d+)[MI]') do -- query length in the alignment
			ql = ql + l
		end
		if bit.band(stack[i][2], 16) == 16 then -- reverse strand
			qb = cigar:match('(%d+)[SH]$')
			if not qb then qb = 0 end
		else -- forward strand
			qb = cigar:match('^(%d+)[SH]')
			if not qb then qb = 0 end
		end
		table.insert(stack[i], tonumber(qb))
		table.insert(stack[i], qb + ql)
		for l in cigar:gmatch('(%d+)[MD]') do -- reference length in the alignment
			tl = tl + l
		end
		table.insert(stack[i], stack[i][4] + tl)
	end
	-- drop overly overlapping hits on the query
	if #stack > 1 then
		local tmp = {stack[1]}
		for i = 2, #stack do
			local dropped = false
			for j = 1, #tmp do
				local b = stack[i][8] > tmp[j][8] and stack[i][8] or tmp[j][8]
				local e = stack[i][9] < tmp[j][9] and stack[i][9] or tmp[j][9]
				if b < e and (e - b) / (stack[i][9] - stack[i][8]) > opt.mask_level then
					dropped = true
					break
				end
			end
			if not dropped then tmp[#tmp+1] = stack[i] end
		end
		stats.n_dropped = stats.n_dropped + (#stack - #tmp)
		stack = tmp
	end
	if opt.is_print then
		for _, v in ipairs(stack) do
			print(v.l)
		end
	end
	for _, v in ipairs(stack) do
		stats.L = stats.L + (v[10] - v[4])
		table.insert(stats.len, v[10] - v[4])
	end
	-- count break points
	if #stack > 1 then
		local b = count_break(stack, opt)
		for i = 1, #b do
			stats.n_b[i] = stats.n_b[i] + (b[i] > 0 and b[i] - 1 or 0)
		end
	else return end
	-- patch small gaps
	if #stack > 1 then
		table.sort(stack, function(a, b) return a[3] < b[3] or (a[3] == b[3] and a[4] < b[4]) end)
		for i = 2, #stack do
			if stack[i][3] == stack[i-1][3] and bit.band(stack[i][2], 16) == bit.band(stack[i-1][2], 16) then -- same chr, same strand
				local gapr = stack[i][4] - stack[i-1][10]
				local gapq = bit.band(stack[i][2], 16) == 0 and stack[i][8] - stack[i-1][9] or stack[i-1][8] - stack[i][9]
				table.insert(stack[i], gapr)
				table.insert(stack[i], gapq)
				if math.abs(gapr) < opt.max_gap and math.abs(gapq) < opt.max_gap then
					stack[i][8] = stack[i][8] < stack[i-1][8] and stack[i][8] or stack[i-1][8]
					stack[i][9] = stack[i][9] > stack[i-1][9] and stack[i][9] or stack[i-1][9]
					stack[i-1][2] = bit.bor(stack[i-1][2], 4)
				end
			end
		end
		local tmp = {}
		for _, v in ipairs(stack) do
			if bit.band(v[2], 4) == 0 then
				tmp[#tmp+1] = v
			end
		end
		stack = tmp
	end
	-- count break points again	
	if #stack > 1 then
		local b = count_break(stack, opt)
		for i = 1, #b do
			stats.n_bg[i] = stats.n_bg[i] + (b[i] > 0 and b[i] - 1 or 0)
		end
	else return end
	-- print the remaining
	--[[
	if not opt.is_print then
		for _, v in ipairs(stack) do
			print(table.concat(v, '\t'))
		end
	end
	]]--
end

function main(arg)
	-- parse the command line
	local opt = {mask_level=0.5, is_print=false, min_q=10, max_gap=500, min_len=200}
	for o, a in os.getopt(arg, 'pq:m:g:l:') do
		if o == 'p' then opt.is_print = true
		elseif o == 'q' then opt.min_q = tonumber(a)
		elseif o == 'm' then opt.mask_level = tonumber(a)
		elseif o == 'g' then opt.max_gap = tonumber(a)
		elseif o == 'l' then opt.min_len = tonumber(a)
		end
	end
	if #arg == 0 then
		print('\nUsage:   luajit ma-stats.lua [options] <bwasw.sam>\n')
		print('Options: -l INT     exclude contigs shorter than INT bp ['..opt.min_len..']')
		print('         -q INT     exclude alignments with maqQ lower than INT ['..opt.min_q..']')
		print('         -m FLOAT   exclude alignments overlapping with a long alignment by FLOAT fraction ['..opt.mask_level..']')
		print('         -g INT     join alignments separated by a gap shorter than INT bp ['..opt.max_gap..']')
		print()
		return
	end
	-- the main part
	local fp, stack, last = io.xopen(arg[1]), {}
	local stats = { L=0, n_un=0, l_un=0, n_dropped=0, len={}, n_b={0,0,0,0,0}, n_bg={0,0,0,0,0}}
	for l in fp:lines() do
		if not l:match('^@') then
			local t = l:split('\t', 10)
			table.remove(t)
			local len = #t[#t] -- sequence length
			if len >= opt.min_len then
				for i = 1, 4 do table.remove(t) end
				t.l = l
				t[#t+1] = len
				t[2], t[4], t[5] = tonumber(t[2]), tonumber(t[4]) - 1, tonumber(t[5]) -- convert to numbers
				if t[1] ~= last then -- change of the query
					analyze(stack, stats, opt)
					stack = {}
					stack[1], last = t, t[1]
				else stack[#stack+1] = t end
			end
		elseif opt.is_print then print(l) end
	end
	analyze(stack, stats, opt)
	fp:close()
	if not opt.is_print then
		-- compute mapped N50
		local l, N50 = 0
		table.sort(stats.len, function(a, b) return b > a end)
		for _, v in ipairs(stats.len) do
			l = l + v
			if l > stats.L / 2 then
				N50 = v
				break
			end
		end
		-- print
		print('CC\tNumber of unmapped contigs: '..stats.n_un)
		print('CC\tTotal length of unmapped contigs: '..stats.l_un)
		print('CC\tNumber of alignments dropped due to excessive overlaps: '..stats.n_dropped)
		print('CC\tMapped contig bases: '..stats.L)
		print('CC\tMapped N50: '..N50)
		print('CC\tNumber of break points: '..stats.n_b[1])
		print('CC\tNumber of Q'..opt.min_q..' break points longer than (0,100,200,500)bp: '
			  .. string.format('(%d,%d,%d,%d)', stats.n_b[2], stats.n_b[3], stats.n_b[4], stats.n_b[5]))
		print('CC\tNumber of break points after closing gaps shorter than '..opt.max_gap..'bp: '..stats.n_bg[1])
		print('CC\tNumber of Q'..opt.min_q..' break points longer than (0,100,200,500)bp after gap closing: '
			  .. string.format('(%d,%d,%d,%d)', stats.n_bg[2], stats.n_bg[3], stats.n_bg[4], stats.n_bg[5]))
	end
end

main(arg)
