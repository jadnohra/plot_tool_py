################################################################################
# Helper: Arguments
################################################################################
import sys

if ( hasattr(sys, 'argv')):
		g_argv = sys.argv
else:
		g_argv = []

g_arg_queried = {}
g_arg_help_on = '-help' in g_argv or '-h' in g_argv or '--h' in g_argv
g_arg_help_arg_count = 0

def print_arg_help():
		global g_arg_help_arg_count
		if (len(g_arg_queried) > g_arg_help_arg_count):
				g_arg_help_arg_count = len(g_arg_queried)
				print 'Help:', ', '.join(g_arg_queried.keys())

def arg_query(keys):
		if (g_arg_help_on and len(keys) and keys[0] not in g_arg_queried):
				g_arg_queried[keys[0]] = ''

def arg_has(keys):
		if (type(keys) is not list):
				keys = [keys]
		arg_query(keys)
		for i in range(len(keys)):
				 if (keys[i] in g_argv):
						return True
		return False
def arg_has_key(keys):
		if (type(keys) is not list):
				keys = [keys]
		arg_query(keys)
		for key in keys:
				ki = g_argv.index(key) if key in g_argv else -1
				if (ki >= 0 and ki+1 < len(g_argv)):
						return True
		return False
def arg_get(keys, dflt):
		if (type(keys) is not list):
				keys = [keys]
		arg_query(keys)
		for key in keys:
				ki = g_argv.index(key) if key in g_argv else -1
				if (ki >= 0 and ki+1 < len(g_argv)):
						return g_argv[ki+1]
		return dflt
def arg_geti(i, dflt):
		arg_query(['index: {}'.format(i)])
		if (i >= len(g_argv)):
				return dflt
		return g_argv[i]

################################################################################
# Helper: Print
################################################################################
import sys

gPrintCol = [ 'default', 'black', 'red', 'green', 'yellow', 'blue', 'magenta', 'cyan', 'white', 'bdefault', 'bblack', 'bred', 'bgreen', 'byellow', 'bblue', 'bmagenta', 'bcyan', 'bwhite'  ]
gPrintColCode = [ "\x1B[0m", "\x1B[30m", "\x1B[31m", "\x1B[32m", "\x1B[33m", "\x1B[34m", "\x1B[35m", "\x1B[36m", "\x1B[37m",
"\x1B[49m", "\x1B[40m", "\x1B[41m", "\x1B[42m", "\x1B[43m", "\x1B[44m", "\x1B[45m", "\x1B[46m", "\x1B[47m", ]
gAltCols = [ gPrintCol.index(x) for x in ['default', 'yellow'] ]

def print_coli(coli):
	coli = coli % len(gPrintCol)
	code = gPrintColCode[coli]
	sys.stdout.write(code)
	#sys.stdout.write('\x1B[{}D'.format(len(code)-3))

def print_col(col):
	print_coli(gPrintCol.index(col))

def printAndChoose(list, postindex = False, forceChoose = False, prefix = ''):
	if (len(list) == 0): return []
	if (len(list) == 1 and forceChoose == False): return list
	for i in range(len(list)):
		print_coli(gAltCols[i % len(gAltCols)])
		if postindex:
			print '{}. {} ({})'.format(prefix, list[i], i+1)
		else:
			print '{}{}. {}'.format(prefix, i+1, list[i])
	print_col('default')
	print '>',
	input_str = raw_input()
	choices = []
	if ('-' in input_str):
		list = input_str.split('-')
		choices = range(int(list[0]), int(list[1])+1)
	elif (',' in input_str):
		choices = [int(x) for x in input_str.split(',')]
	else:
		if len(input_str):
			choices.append(int(input_str))
	choices = [i-1 for i in choices]
	return choices

################################################################################
# Helper: Miscellaneous
################################################################################
def selRange(input_str):
	def handleRange(rng):
		return range(int(rng.split(':')[0].split('-')[0]), int(rng.split(':')[0].split('-')[1])+1, int(rng.split(':')[1]) if ':' in rng else 1)
	return sum( [handleRange(x) if ('-' in x) else [int(x)] for x in input_str.split(',')], [])

################################################################################
# Helper: Math
################################################################################
import numpy
import math

def rand_between(a, b):
		return a+(numpy.random.rand() * (b-a))
def rand_ampl(ampl):
		return 2.0*(numpy.random.rand()-0.5) * ampl
def rand_sgn():
		return m_sgn(rand_ampl(1.0))
def rand_ind(ln):
		return int(numpy.clip(round(numpy.random.rand()*float(ln)), 0, ln-1))

def vec_prec_str(v, prec):
		if (prec <= 0) or (prec is None):
				return v
		fmt = '{{:.{}g}}'.format(prec); return ','.join([fmt.format(x) for x in v]);
def vec_prec(v, prec):
		return [float(x) for x in vec_prec_str(v, prec).split(',')]

def rand_ampl_prec(ampl, prec=2):
		return vec_prec([rand_ampl(ampl)], prec)[0]

def is_float(s):
		try:
				float(s)
				return True
		except ValueError:
				return False

################################################################################
# plot_tool
################################################################################

from sympy import *
import numpy
import matplotlib.pyplot as pyplot
import matplotlib.cm
import mpl_toolkits.mplot3d.axes3d as axes3d
import matplotlib.animation as animation
import re
import time
import threading, multiprocessing
import os
import traceback
import fcntl, os, sys
import scipy.cluster.vq

g_ctx = {
	'func_str' : None, 'func_lambd': None,
	'symb_consts' : {},
	'val_consts' : {},
	'symb_vars' : {},
	'ranges' : {},
	'range_offsets' : {},
	'parse_func_re' : re.compile(ur'((?:\d+\.*(?:0?e-\d+)?\d*(?!\w)|\w+)\b)(?!\s*\()'),
	'parse_number_re' : re.compile(ur'\d+\.*(?:0?e-\d+)?\d*(?!\w)'),
	'parse_symb_re' : re.compile(ur'\w+'),
	'parse_arithm_re' : re.compile(ur'(\+|\-|\*|\/|\\|\^|\.)'),
	'parse_constants' : 'pi'.split(','),
	'known_vars' : 'x,y,z,w,t,p,q,r,s,t,u,v'.split(','),
	'known_consts' : 'a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p'.split(','),
	'range_scl' : 1.0, 'resolution_scl' : 1.0,
	'style' : 'surface',
	'cmd' : '', 'last_cmd' : '', 'replot' : False,
	'anim_t' : [0.0], 'anim_t_dir' : [1.0], 'anim_f_lambda' : None,
	'anim_vrs': None, 'anim_objs' : None, 'anim_axs': None, 'anim_t_speed' : 1.0,
	'anim_synced' : True, 'anim_t_n':0, 'anim_paused' : False, 'anim_one_step' : False,
	'family_rms' : 1.0, 'countour_count': 10,
	'case_swap' : False,
	'zlim_mode' : 'za', 'zlim_scl' : 1.0,
	'rand_radius' : 1.0, 'rand_digits' : 2,
	'mouse_coords' : None, 'cam_coords' : None,
	'cc_str' : None, 'scatter':None, 'init_scatter':False, 'show_scatter':True,
	'line_width':0.75,
	'script':None, 'init_script':False,
}
def get_ctx():
	global g_ctx; return g_ctx;
def as_list(names):
	return names.split(',') if (type(names) is not list) else names
def has_symb(name):
	ctx = get_ctx(); return name in ctx['symb_consts'] or name in ctx['symb_vars'];
def has_const(name):
	ctx = get_ctx(); return name in ctx['symb_consts'];
def def_const(names, val = 0.5):
	ctx = get_ctx()
	names = as_list(names)
	for name in names:
		for d in (ctx['symb_consts'],ctx['symb_vars']):
			d.pop(name, None)
		ctx['symb_consts'][name] = Symbol(name)
		ctx['val_consts'][name] = val
def val_of_const(name, def_val = 0.5):
	ctx = get_ctx()
	if not has_const(name):
		def_const(name, def_val)
	return ctx['val_consts'][name]
def set_const(name, val = 0.5):
	ctx = get_ctx()
	if not has_const(name):
		def_const(name, val)
	else:
		ctx['val_consts'][name] = val
def def_var(names):
	ctx = get_ctx()
	names = as_list(names)
	for name in names:
		for d in (ctx['symb_consts'],ctx['symb_vars']):
			d.pop(name, None)
		ctx['symb_vars'][name] = Symbol(name)
def def_range(names, rng):
	ctx = get_ctx()
	names = as_list(names)
	for name in names:
		ctx['ranges'][name] = [float(x) for x in as_list(rng)]
def has_range(name):
	ctx = get_ctx(); return name in ctx['ranges'];
def range_of(name, scl=1.0, res=1.0):
	ctx = get_ctx(); rng = ctx['ranges'][name]
	rad = (rng[1]-rng[0])/2.0; o = ctx['range_offsets'].get(name, 0.0);
	return [o-(rad*scl), o+(rad*scl), rng[2]*res]
def parse_func_str(func_str):
	# e.g: 'cos(xy+sin (x_z^2))-0.2*h+1.0e-10+2+g2+pi*5'
	ctx = get_ctx()
	symbols = []; numbers = []; constants = [];
	matches = re.findall(ctx['parse_func_re'], func_str)
	for m in matches:
		#print m
		subm = re.findall(ctx['parse_number_re'], m)
		if len(subm) == 1 and subm[0] == m:
			numbers.append(m)
			continue
		subm = re.findall(ctx['parse_symb_re'], m)
		if len(subm) == 1 and subm[0] == m:
			if (m.lower() in ctx['parse_constants']):
				constants.append(m)
			else:
				symbols.append(m)
			continue
	return (symbols, numbers, constants)
def form_func_str(func_str, cprefix = 'C'):
	(symbols, numbers, constants) = parse_func_str(func_str)
	form_func = func_str; form_consnts = [];
	for num in numbers:
		if ( floor(float(num)) != float(num)):
			form_c = '{}{}'.format(cprefix, len(form_consnts)); form_consnts.append(form_c);
			form_func = form_func.replace(num, form_c)
	return form_func, form_consnts
def eval_symb_str(str, symbs, do_expand):
	for symb in symbs:
		exec('{} = Symbol("{}")'.format(symb, symb))
	ev = eval(str); return expand(ev) if do_expand else ev;
def get_range_scale():
	ctx = get_ctx(); return ctx['range_scl'];
def scale_range(scl):
	ctx = get_ctx(); ctx['range_scl'] = scl;
def get_resolution_scale():
	ctx = get_ctx(); return ctx['resolution_scl'];
def scale_resolution(scl):
	ctx = get_ctx(); ctx['resolution_scl'] = scl;
def get_anim_t_speed():
	ctx = get_ctx(); return ctx['anim_t_speed'];
def set_anim_t_speed(val):
	ctx = get_ctx(); ctx['anim_t_speed'] = val;
def set_style(style):
	ctx = get_ctx(); ctx['style'] = style;
def get_func_val(coords):
	return get_ctx()['func_lambd'](*coords)
def get_func_val_str(coords):
	return 'f({},{})= {}'.format( *(coords[:2] + [get_func_val(coords[:2])]) )
def to_scatter(data):
	data = data if data is not None else []
	return ([x[0] for x in data], [x[1] for x in data], [(x[2] if len(x)>2 else None) for x in data])
def _at(coords):
	return 'f({},{})= {}'.format( *(coords[:2] + [get_func_val(coords[:2])]) )
def _p(pi):
	return [get_ctx()['scatter'][x][pi] for x in range(3)]
def __around(coords, rad, cnt):
	def _around_func(coords, rad, cnt, data):
		_data = []; rad = [-rad, rad] if (type(rad) is not list) else rad;
		for i in range(cnt):
			dc = [rand_between(*rad) for x in range(2)]; pt = [coords[x]+dc[x] for x in range(2)];
			_data.append([pt[0], pt[1], get_func_val(pt)])
		data.extend(sorted(_data, key = lambda el:el[2]))
	nthreads = max(1, multiprocessing.cpu_count()-1); threads = []; datas = [[] for x in range(nthreads)]
	for ti in range(nthreads):
		tn = max(1, int(round( float(cnt) /float(nthreads))))
		td = threading.Thread(target=_around_func, args=(coords, rad, tn, datas[ti]))
		threads.append(td); td.setDaemon(True); td.start();
	for td in threads:
		td.join()
	all_data = []
	for data in datas:
		all_data.extend(data)
	all_data = sorted(all_data, key = lambda el:el[2])
	return all_data
def _around(coords, rad, cnt):
	t0 = time.clock()
	data = __around(coords, rad, cnt)
	print '({:.2f} secs.)'.format(time.clock() - t0);
	return data
def _form(key='func_sym_exp'):
	return to_std_exp(form_func_str(str(get_ctx()[key]))[0])
def handle_console_input_init():
	fd = sys.stdin.fileno()
	fl = fcntl.fcntl(fd, fcntl.F_GETFL)
	fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)
def handle_console_input():
	try:
		input = sys.stdin.readline()
		#print 'echo:', input,
		return input
	except:
		return ''
def process_cmd(cmd_str, full = True, allow_close_plot = True):
	def cmd_close_plot(replot):
		if (allow_close_plot):
			close_plot(replot)
	def is_float(s):
		try:
			float(s)
			return True
		except ValueError:
			return False
	ctx = get_ctx()
	repeat = True; pre_mul = 1.0; reg_cmd_str = cmd_str; cmd_print = True; ri = 0;
	try:
		while(repeat and (len(cmd_str) or full)):
			splt = cmd_str.split(' '); cmd = splt[0];
			if len(cmd_str) and cmd_print and arg_has('-print_cmd'):
				print 'Command: [{}]'.format(cmd_str)
			repeat = False
			if (cmd == 'C'):
				cmd_str = handle_console_input().strip(); full = cmd_str != ''; repeat = True;
			elif (cmd == 'zo'):
				scale_range(get_range_scale()*2.0*pre_mul); scale_resolution(get_resolution_scale()*2.0*pre_mul);
				ctx['family_rms'] = ctx['family_rms']*2.0*pre_mul; cmd_close_plot(True);
			elif (cmd == 'zi'):
				scale_range(get_range_scale()/(2.0*pre_mul)); scale_resolution(get_resolution_scale()/(2.0*pre_mul));
				ctx['family_rms'] = ctx['family_rms']/(2.0*pre_mul); cmd_close_plot(True);
			elif (cmd == 'z' and full):
				fact = float(splt[1])
				scale_range(get_range_scale()/(fact*pre_mul)); scale_resolution(get_resolution_scale()/(fact*pre_mul));
				ctx['family_rms'] = ctx['family_rms']/(2.0*pre_mul); cmd_close_plot(True);
			elif (cmd == 'zr'):
				scale_range(1.0); scale_resolution(1.0); ctx['family_rms'] = 1.0; cmd_close_plot(True);
			elif (cmd == 'ro'):
				scale_range(get_range_scale()*2.0*pre_mul); cmd_close_plot(True);
			elif (cmd == 'ri'):
				scale_range(get_range_scale()/(2.0*pre_mul)); cmd_close_plot(True);
			elif (cmd == 'co'):
				ctx['countour_count'] = max(2, ctx['countour_count']/2); cmd_close_plot(True);
			elif (cmd == 'ci'):
				ctx['countour_count'] = ctx['countour_count']*2; cmd_close_plot(True);
			elif (cmd == 'reso'):
				scale_resolution(get_resolution_scale()*2.0*pre_mul); cmd_close_plot(True);
			elif (cmd == 'resi'):
				scale_resolution(get_resolution_scale()/(2.0*pre_mul)); cmd_close_plot(True);
			elif (cmd == 'rmso'):
				ctx['family_rms'] = ctx['family_rms']*2.0*pre_mul; cmd_close_plot(True);
			elif (cmd == 'rmsi'):
				ctx['family_rms'] = ctx['family_rms']/(2.0*pre_mul); cmd_close_plot(True);
			elif (cmd == 'w'):
				set_style('wire'); cmd_close_plot(True);
			elif (cmd == 's'):
				set_style('surface'); cmd_close_plot(True);
			elif (cmd == 'c' and full):
				set_style('contour'); cmd_close_plot(True);
			elif (cmd == 'x'):
				set_style('xray'); cmd_close_plot(True);
			elif (cmd == 'zl'):
				ctx['zlim_mode'] = 'zl'; cmd_close_plot(True);
			elif (cmd == 'zu'):
				ctx['zlim_mode'] = 'zu'; cmd_close_plot(True);
			elif (cmd == 'zm'):
				ctx['zlim_mode'] = 'zm'; cmd_close_plot(True);
			elif (cmd == 'zx'):
				ctx['zlim_mode'] = 'zx'; cmd_close_plot(True);
			elif (cmd == 'zzi'):
				ctx['zlim_scl'] = ctx['zlim_scl']*0.5; cmd_close_plot(True);
			elif (cmd == 'zzo'):
				ctx['zlim_scl'] = ctx['zlim_scl']*2.0; cmd_close_plot(True);
			elif (cmd == 'zzr'):
				ctx['zlim_scl'] = 1.0; cmd_close_plot(True);
			elif (cmd == 'f'):
				set_style('family'); cmd_close_plot(True);
			elif (cmd == 'lw' and full):
				if (float(splt[1]) > 0):
					ctx['line_width'] = float(splt[1]); cmd_close_plot(True);
			elif (cmd == 'lwo'):
				ctx['line_width'] = ctx['line_width']*1.5; cmd_close_plot(True);
			elif (cmd == 'lwi'):
				ctx['line_width'] = ctx['line_width']/1.5; cmd_close_plot(True);
			elif (cmd == 'y'):
				ctx['anim_synced'] = not (ctx['anim_synced'])
			elif (cmd == 'Z'):
				if ('Z' in ctx):
					Z = ctx['Z']; Zu = {};
					for row in Z:
						for z in row:
							if (z in Zu):
								Zu[z] = Zu[z]+1
							else:
								Zu[z] = 1
					pi = 0
					print '----'
					for k in sorted(Zu.keys(), cmp=lambda k1,k2: Zu[k2]-Zu[k1] if (Zu[k1] != Zu[k2]) else (-1 if k1<k2 else (1 if k1>k2 else 0))):
						print ' {} x {}{}'.format(Zu[k], k, '' if (pi % 6) else '\n'),
						pi = pi+1
					print ''
			elif (cmd == '>'):
				set_anim_t_speed(get_anim_t_speed()*2.0)
				if ctx['style'] in 'f,fam,family,families'.split(','):
					cmd_close_plot(True);
			elif (cmd == '<'):
				set_anim_t_speed(get_anim_t_speed()*0.5)
				if ctx['style'] in 'f,fam,family,families'.split(','):
					cmd_close_plot(True);
			elif (cmd == '[?'):
				print '------'
				for v in ctx['func_vrs']:
					print '{}: {}'.format(v, range_of(v, get_range_scale(), get_resolution_scale()))
			elif (cmd.strip() == '' and full):
				cmd_str = ctx['last_cmd']; splt = cmd_str.split(' '); cmd = splt[0]; repeat = True; reg_cmd_str = cmd_str;
			elif (cmd_str == ' ' and (not full) and ri == 0):
				ctx['anim_one_step'] = True
			elif (cmd == 'p'):
				ctx['anim_paused'] = not ctx['anim_paused']
			elif (cmd == 'P'):
				ctx['show_scatter'] = not ctx['show_scatter']; cmd_close_plot(True);
			elif (cmd == 'S'):
				ctx['case_swap'] = not ctx['case_swap'];
				ctx['symb_consts'] = {}; ctx['val_consts'] = {}; ctx['symb_vars'] = {}; #ctx['ranges'] = {};
				cmd_close_plot(True);
			elif (cmd.startswith(':?')):
				print '------'
				for k,v in ctx['val_consts'].items():
					 print '{}: {}'.format(k,v)
			elif (cmd == 'R'):
				ampl = ctx['rand_radius']; prec = ctx['rand_digits'];
				print 'Rand: '
				for k in ctx['val_consts'].keys():
					v = ctx['val_consts'][k] = vec_prec([ 2.0*(numpy.random.rand()-0.5) * ampl], prec)[0]
					print '{}: {}, '.format(k,v),
				print ''
				cmd_close_plot(True);
			elif (cmd == 'ra' and full):
				ctx['rand_radius'] = float(splt[1])
			elif (cmd == 'rd' and full):
				ctx['rand_digits'] = int(splt[1])
			elif (cmd == 'rr' and full):
				for k in ctx['val_consts'].keys():
					ctx['val_consts'][k] = 0.5
			elif (full and cmd == ':*'):
				vals = splt[1:]
				keys = sorted(ctx['symb_consts'])
				for i in range(len(vals)):
					if is_float(vals[i]):
						print '{} : {}'.format(keys[i], float(vals[i]));
						set_const(keys[i], vals[i]);
				cmd_close_plot(True);
			elif (full and cmd.startswith(':')):
				if (cmd_str.count(':') == 1):
					name = cmd[1:]; val = float(splt[1]); set_const(name, val); print cmd_str; cmd_close_plot(True);
				else:
					print 'Ignored: [{}]'.format(cmd_str)
			elif (full and cmd.startswith('o')):
				offs = ctx['range_offsets']
				if (cmd == 'o'):
					print '------'
					for k,v in  offs.items():
						 print '{}: {}'.format(k,v)
				elif len(splt)>1:
					names = cmd[1:].split(',')
					op = splt[1][0]
					if op in '+,*,/'.split(','):
						val = float(splt[1][1:])
						for name in names:
							offs[name] = eval('{} {} {}'.format(offs.get(name, 0.0), op, val))
					else:
						for name in names:
							offs[name] = float(splt[1])
					cmd_close_plot(True);
			elif (cmd in ['cx', 'cy', 'cc']) or (full and cmd.startswith('cp')):
				if full:
					pi = int(cmd[2:]); coords = [ctx['scatter'][i][pi] for i in range(3)]; vis = [0,1];
				else:
					vis = {'cx':[0], 'cy':[1], 'cc':[0,1] }[cmd]
					coords = ctx['mouse_coords'][:2]
				ctx['last_cc'] = coords
				for vi in vis:
					name = ctx['func_vrs'][vi]; off = coords[vi];
					ctx['range_offsets'][name] = off; print '{}: {}'.format(name, off);
				cmd_close_plot(True);
			elif (cmd == 'D'):
				center = [0.0, 0.0]
				for vi in range(2):
					name = ctx['func_vrs'][vi]; center[vi] = ctx['range_offsets'].get(name, 0.0)
				rng = range_of(ctx['func_vrs'][0]); rad = (rng[1]-rng[0])/1.0;
				lowest = __around(center, rad, 1.e5)[0]
				print _at(center); print _at(lowest);
				if (get_func_val(lowest[:2]) < get_func_val(center)):
					for vi in range(2):
						name = ctx['func_vrs'][vi]; off = lowest[vi];
						ctx['range_offsets'][name] = off; print '{}: {}'.format(name, off);
						cmd_close_plot(True);
			elif (cmd == '='):
				print get_func_val_str(ctx['mouse_coords'])
			elif (full and cmd.startswith('@p')):
				pi = int(cmd[2:]); coords = [ctx['scatter'][i][pi] for i in range(3)];
				print get_func_val_str(coords)
			elif (full and cmd.startswith('@')):
				coords = [float(x) for x in cmd[1:].split(',')]; print get_func_val_str(coords);
			elif (cmd == '@cc'):
				print get_func_val_str(ctx['last_cc'])
			elif (is_float(cmd) and len(splt) > 1):
				pre_mul = float(cmd); splt.pop(0);
				if (len(splt) and len(splt[0])):
					cmd_str = ' '.join(splt); repeat = True; cmd_print = False;
				else:
					return False
			else:
				return False
			ri = ri+1
		ctx['cmd'] = ''
		if (reg_cmd_str is not None):
			if arg_has('-print_cmd'):
				print 'Executed: [{}]'.format(reg_cmd_str)
			ctx['last_cmd'] = reg_cmd_str
	except:
		traceback.print_exc()
		ctx['cmd'] = ''; ctx['last_cmd'] = '';
	return True
def plot_on_key(event):
	ctx = get_ctx()
	if (event.key == 'enter'):
		process_cmd(ctx['cmd'])
		ctx['cmd'] = ''
	else:
		key = event.key
		ctx['cmd'] = ctx['cmd'] + reduce(lambda k,(x,y):k.replace(x, y), [key, ('ctrl+',''),('cmd+',''),('alt+',''),('space',' ')])
		if (not process_cmd(ctx['cmd'], False)) and len(ctx['cmd']) > 2:
			print ctx['cmd']
def plot_on_motion(event):
	try:
		ctx = get_ctx();
		ctx['mouse_coords'] = [float(x.split('=')[1]) for x in pyplot.gca().format_coord(event.xdata, event.ydata).split(',')]
	except:
		return
def plot_on_release(event):
	ctx = get_ctx();
	ctx['cam_coords'] = (pyplot.gca().azim, pyplot.gca().elev)
def close_plot(replot):
	get_ctx()['replot'] = replot; get_ctx()['cmd'] = ''; pyplot.close();
def limit_z(Z, Xl, Yl, Zr):
	ctx = get_ctx(); mode = ctx['zlim_mode']; zscl = ctx['zlim_scl'];
	zl = None; cl = None;
	if (mode == 'zl'):
		zl = [Zr[0], Zr[0]+zscl*(Xl[1]-Xl[0])]
	elif (mode == 'zhl'):
		zl = [Zr[0], Zr[0]+zscl*(Xl[1]-Xl[0])/2.0]
	elif (mode == 'zu'):
		zl = [Zr[1]-zscl*(Xl[1]-Xl[0]), Zr[1]]
	elif (mode == 'zhl'):
		zl = [Zr[1]-zscl*(Xl[1]-Xl[0])/2.0, Zr[1]]
	#print zl
	if (zl is not None) and (cl is None):
		zll = zl[1]-zl[0]; cl = [zl[0]-0.01*zll, zl[1]+0.01*zll]
	if (cl is not None) and (Zr[0] < cl[0] or Zr[1] > cl[1]):
		for i in range(len(Z)):
			Z[i] = numpy.clip(Z[i], cl[0], cl[1])
	return zl
def range_put(r, v):
	r[0] = min(r[0], v); r[1] = max(r[1], v);
def anim_calc_data(func_lambd, vrs, vi, ts):
	ctx = get_ctx(); rng_scl = ctx['range_scl']; res_scl = ctx['resolution_scl'];
	X = numpy.arange(*range_of(vrs[vi], rng_scl, res_scl))
	Z = []; Zr = [float('inf'), float('-inf')];
	vals = [0.0]*(len(vrs))
	ti = 0
	init_val_rng = range(len(vrs)) if (len(ts) == len(vrs)) else [x for x in range(len(vrs)) if (x!=vi)]
	for i in init_val_rng:
		yrng = [rng_scl*x for x in range_of(vrs[i])]; vals[i] = yrng[0] + ((1.0)*ts[ti] * (yrng[1]-yrng[0])); ti=ti+1;
	for i in range(len(X)):
		vals[vi] = X[i]; z = func_lambd(*vals); Z.append( z ); range_put(Zr, z);
	return X,Z,vals
def anim_plot_update(num):
	ctx = get_ctx(); pobjs = ctx['anim_objs']; axs = ctx['anim_axs']; vrs = ctx['anim_vrs'];
	if ctx['anim_paused'] and (not ctx['anim_one_step']):
		return
	ctx['anim_one_step'] = False
	for vi in range(len(vrs)):
		X,Z,vals = anim_calc_data(ctx['anim_f_lambda'], vrs, vi, ctx['anim_t'])
		pobjs[vi].set_data(X,Z)
		ax = axs[vi]
		ax.relim(); ax.autoscale_view(True, False, True);
		locs = ax.get_xticks()[:-(len(vrs)-1)]; nticks = [x for x in locs] + [vals[i] for i in range(len(vals)) if i != vi]; ax.set_xticks(nticks);
		labels = ax.get_xticklabels();
		ax.set_xticklabels(nticks[:-(len(vrs)-1)] + [x for x in vrs if x!=vrs[vi]]);
	for ti in range(ctx['anim_t_n']):
		t = ctx['anim_t'][ti] + 0.05*ctx['anim_t_dir'][ti]*ctx['anim_t_speed']
		if (t > 1.0):
			t = 1.0; ctx['anim_t_dir'][ti] = -1.0
		if (t < 0.0):
			t = 0.0; ctx['anim_t_dir'][ti] = 1.0
			ctx['anim_t'][ti] = t
		else:
			ctx['anim_t'][ti] = t
			break
def to_std_delims(frag):
	return frag.replace('[','(').replace(']',')')
def to_sympy_exp(frag):
	return frag.replace('^', '**')
def to_std_exp(frag):
	return frag.replace('**', '^')
def format_func_str(func_str):
	return to_sympy_exp(to_std_delims(func_str))
def plot(_func_str, rg = None, rg_x = None, rg_y =None, res = None, res_x = None, res_y = None, _rg = None, _rg_x = None, _rg_y = None):
	ctx = get_ctx()
	func_str = format_func_str(_func_str)
	(func_symbols, func_numbers, func_constants) = parse_func_str(func_str)
	func_symbols = sorted(func_symbols)
	vrs = []; constnts = []; func_lambd = None; pre_func_sym = None; func_sym = None;
	for symb in func_symbols:
		is_var = (symb.upper() != symb)
		if ( is_var if (not ctx['case_swap']) else (not is_var) ):
			if (not has_symb(symb)):
				print 'Registered variable', symb
				def_var(symb)
				if (not has_range(symb)):
					def_range(symb, [-1.0,1.0,0.1])
			if (symb not in vrs):
				vrs.append(symb)
		else:
			if (not has_symb(symb)):
				print 'Registered constant', symb
				def_const(symb)
			if (symb not in constnts):
				constnts.append(symb)
	if ctx['init_script']:
		ctx['init_script'] = False
		print 'Executing script...'
		for il in range(len(ctx['script'])):
			print ' {}. {}'.format(il+1, ctx['script'][il])
			process_cmd(ctx['script'][il], True, False)
		print 'done'
	func = eval_symb_str(func_str, vrs + constnts, False)
	pre_func_sym = func
	func_sub = func
	for cnt in constnts:
		func_sub = func_sub.subs(cnt, val_of_const(cnt))
	func_sym = func_sub
	func_lambd = lambdify(vrs, expand(func_sub), "numpy")
	rng_scl = ctx['range_scl']; res_scl = ctx['resolution_scl'];
	fig = None; ctx['Z'] = None;
	fallback_style = ctx['style']
	if fallback_style in 'w,wire,wireframe,s,surf,surface,contour'.split(','):
		if (len(vrs) > 2):
			fallback_style = 'x'
		if (len(vrs) == 1):
			fallback_style = '2d'
	if fallback_style in '2d'.split(','):
		if ctx['cc_str'] is not None:
			cc = [float(x) for x in ctx['cc_str'].split(',')]
			ctx['last_cc'] = cc; ctx['cc_str'] = None;
			for vi in range(len(vrs)):
				ctx['range_offsets'][vrs[vi]] = cc[vi]
		Xr = range_of(vrs[0], rng_scl, res_scl);
		X = numpy.arange(*Xr)
		Z = []; Zr = [float('inf'), float('-inf')];
		for i in range(len(X)):
			row = []
			z = func_lambd(X[i]); row.append( z ); range_put(Zr, z);
			Z.append(row)
		zlim = limit_z(Z, Xr, Xr, Zr)
		fig, ax = pyplot.subplots(len(vrs));
		if (zlim is not None):
			ax.set_ylim(zlim)
		ax.plot(X, Z);
	if fallback_style in 'w,wire,wireframe,s,surf,surface,contour'.split(','):
		if ctx['cc_str'] is not None:
			cc = [float(x) for x in ctx['cc_str'].split(',')]
			ctx['last_cc'] = cc; ctx['cc_str'] = None;
			for vi in range(len(vrs)):
				ctx['range_offsets'][vrs[vi]] = cc[vi]
		Xr = range_of(vrs[0], rng_scl, res_scl); Yr = range_of(vrs[1], rng_scl, res_scl);
		X, Y = numpy.meshgrid(numpy.arange(*Xr), numpy.arange(*Yr))
		Z = []; Zr = [float('inf'), float('-inf')];
		for i in range(len(X)):
			row = []
			for j in range(len(X[0])):
				z = func_lambd(X[i,j], Y[i,j]); row.append( z ); range_put(Zr, z);
			Z.append(row)
		if fallback_style != 'contour':
			fig = pyplot.figure()
			zlim = limit_z(Z, Xr, Yr, Zr)
			ax = axes3d.Axes3D(fig)
			if (ctx['cam_coords'] is not None):
				ax.view_init(ctx['cam_coords'][1], ctx['cam_coords'][0])
			if (zlim is not None):
				ax.set_zlim3d(zlim)
			if ctx['style'] in 'w,wire,wireframe'.split(','):
				ax.plot_wireframe(X,Y,Z,cmap = matplotlib.cm.coolwarm, linewidth=ctx['line_width'], color='k' if arg_has('-paper') else 'b')
			else:
				ax.plot_surface(X,Y,Z,cmap = matplotlib.cm.coolwarm, rstride = 1, cstride = 1, alpha=1.0, linewidth=ctx['line_width'])
			if ctx['show_scatter'] and  (ctx['scatter'] is not None):
				sX,sY,sZ = (ctx['scatter'])
				if (ctx['init_scatter']):
					ctx['init_scatter'] = False
					for i in range(len(sZ)):
						sZ[i] = sZ[i] if (sZ[i] is not None) else func_lambd(sX[i], sY[i])
				ax.scatter(sX,sY,sZ, s=32, picker=True, c='k' if arg_has('-paper') else 'b', facecolors='none', alpha=1.0)
			ctx['Z'] = Z
			fig.canvas.mpl_connect('button_release_event', plot_on_release)
		else:
			levels = {};
			for i in range(len(X)):
				for j in range(len(X[0])):
					x,y,z = X[i,j], Y[i,j], Z[i][j]
					levels[z] = z
			centroids = scipy.cluster.vq.kmeans(numpy.asarray(levels.keys()), ctx['countour_count'])[0]
			centroids = sorted(centroids)
			contours = {}
			for ci in range(len(centroids)):
				contours[centroids[ci]] = ([],[], ci)
			for i in range(len(X)):
				for j in range(len(X[0])):
					x,y,z = X[i,j], Y[i,j], Z[i][j]
					closest_ci = 0; closest_dist = float('inf');
					#print '{},{},{}'.format(x,y,z)
					for ci in range(len(centroids)):
						cent = centroids[ci]
						if (abs(z-cent) < closest_dist):
							#print '{}-{}'.format(abs(z-cent), ci)
							closest_dist = abs(z-cent); closest_ci = ci;
					#if (closest_ci == 0):
					#	print '{},{},{}'.format(x,y,z)
					#print '{} -> {} {}'.format(z, closest_ci, centroids[closest_ci])
					cont = contours[centroids[closest_ci]];	cont[0].append(x); cont[1].append(y);
			fig, ax = pyplot.subplots(1)
			ax.set_xlim(Xr[0], Xr[1]); ax.set_ylim(Yr[0], Yr[1]);
			for cont in [x for x in contours.values() if len(x[0])]:
				shade = float(cont[2])/float(len(centroids)); sp = [0.1*x for x in numpy.random.rand(3)];
				ax.scatter(cont[0], cont[1], c=(min(shade+sp[0],1.0),min(sp[2],1.0),1.0-min(shade+sp[1],1.0)))
				#break
	elif fallback_style in 'f,fam,family,families'.split(','):
		def rms(x):
			return numpy.sqrt(x.dot(x)/x.size)
		def family_curve_check(func_lambd, vrs, vi, anim_ts, families, fam_rms):
			X,Z,vals = anim_calc_data(func_lambd, vrs, vi, anim_ts)
			has_curve = False; Za = numpy.asarray(Z);
			for curve in families[vi]:
				if rms(Za-curve) < fam_rms:
					has_curve = True; break;
			if not has_curve:
				families[vi].append(Za)
		def fam_step_times(anim_ts):
			could_step = False
			for ti in range(len(anim_ts)):
				anim_ts[ti] = anim_ts[ti] + 0.05*get_anim_t_speed()
				if (anim_ts[ti] <= 1.0):
					could_step = True; break;
				else:
					anim_ts[ti] = 0.0
			return could_step
		def anim_family_curve_check(func_lambd, vrs, vi, families, fam_rms):
			could_step = True; anim_ts = [0.0]*(len(vrs)-1);
			while could_step:
				family_curve_check(func_lambd, vrs, vi, anim_ts, families, fam_rms)
				could_step = fam_step_times(anim_ts)
		families = [[] for x in vrs];
		use_threads = not arg_has('-no_threads');
		print 'Finding families...',; sys.stdout.flush(); t0 = time.clock();
		if (use_threads):
			threads = []
			for vi in range(len(vrs)):
				td = threading.Thread(target=anim_family_curve_check, args=(func_lambd, vrs, vi, families, ctx['family_rms']))
				threads.append(td); td.setDaemon(True); td.start();
			for td in threads:
				td.join()
		else:
			could_step = True; anim_ts = [0.0]*(len(vrs));
			while could_step:
				for vi in range(len(vrs)):
					family_curve_check(func_lambd, vrs, vi, anim_ts, families, ctx['family_rms'])
				could_step = fam_step_times(anim_ts)
		print ' ({:.2f} secs.)'.format(time.clock() - t0); sys.stdout.flush();
		fig, axi = pyplot.subplots(len(vrs))
		for vi in range(len(vrs)):
			axi[vi].set_xlabel(vrs[vi]); axi[vi].xaxis.set_label_coords(1.05, -0.025);
			for Z in [x.tolist() for x in families[vi]]:
				X = numpy.arange(*range_of(vrs[vi], rng_scl, res_scl))
				axi[vi].plot(X, Z)
	elif fallback_style in 'x,xray'.split(','):
		ctx['anim_f_lambda'] = func_lambd; ctx['anim_vrs'] = vrs;
		ctx['anim_t_n'] = (len(vrs)) if ctx['anim_synced'] else (len(vrs)-1);
		ctx['anim_t'] = [0.0]*(ctx['anim_t_n']); ctx['anim_t_dir'] = [1.0]*(ctx['anim_t_n']);
		ctx['anim_objs'] = []; ctx['anim_axs'] = [];
		fig, axi = pyplot.subplots(len(vrs));
		if (type(axi) is not list):
			axi = [axi]
		for vi in range(len(vrs)):
			X,Z,vals = anim_calc_data(func_lambd, vrs, vi, ctx['anim_t'])
			ax = axi[vi]; pobj = ax.plot(X, Z); ax.set_xlabel(vrs[vi]); ax.xaxis.set_label_coords(1.05, -0.025);
			ctx['anim_objs'].append(pobj[0]); ctx['anim_axs'].append(ax);
			locs = ax.get_xticks(); nticks = [x for x in locs] + ([0.0]*(len(vrs)-1)); ax.set_xticks(nticks);
			labels = ax.get_xticklabels();
			ax.set_xticklabels(nticks[:-(len(vrs)-1)] + [x for x in vrs if x!=vrs[vi]]);
		anim = animation.FuncAnimation(fig, anim_plot_update, None, None, interval=30, blit=False)
	ctx['func_str'] = _func_str; ctx['func_vrs'] = vrs;
	ctx['func_lambd'] = func_lambd; ctx['pre_func_sym'] = pre_func_sym; ctx['func_sym'] = func_sym;
	ctx['pre_func_sym_exp'] = expand(pre_func_sym) if pre_func_sym else None
	ctx['func_sym_exp'] = expand(func_sym) if func_sym else None
	if (fig):
		fig.canvas.mpl_connect('key_press_event', plot_on_key)
		fig.canvas.mpl_connect('motion_notify_event', plot_on_motion)
		pyplot.grid()
		pyplot.show()
def replot():
	ctx = get_ctx()
	if 'func_str' in ctx:
		plot(ctx['func_str'])
def init():
	for k in pyplot.rcParams.keys():
		if k.startswith('keymap.'):
			pyplot.rcParams[k] = ''
	pyplot.rcParams['keymap.quit'] = 'q,escape'
def post_init():
	handle_console_input_init
def def_repo_fp():
	return 'plot_tool_repo.txt'
def read_repo(fp = def_repo_fp(), desc_filter = [], with_scatter = False):
	repo = []
	funci = 0
	with open(fp, "r") as ifile:
		desc = ''; multiline = []; multiscatter = []; multiscript = [];
		for line in ifile.readlines():
			line = line.strip();
			if len(line):
				if line.startswith('###'):
					multiscript.append(line[3:].strip())
				elif line.startswith('##'):
					multiscatter.append(line[2:].strip())
				elif line.startswith('#'):
					if len(multiline) == 0:
						desc = '{}'.format(line[1:].strip());
				else:
					if line.endswith('..'):
						multiline.append(line[:-2])
					else:
						line = ''.join(multiline + [line])
						if all([(pat.lower() in desc.lower()) for pat in desc_filter]) and ((not with_scatter) or (len(multiscatter))):
							scatter = eval('[{}]'.format(''.join(multiscatter))) if len(multiscatter) else None
							repo.append( { 'func':line, 'desc':desc, 'i':funci, 'scatter':scatter, 'script':multiscript} )
							funci = funci+1
						desc = ''; multiline = []; multiscatter = []; multiscript = [];
	return repo
def add_to_repo(func_str, func_desc, fp = def_repo_fp()):
	if (not os.path.isfile(fp)):
		with open(fp, "w") as ofile:
			ofile.write('#\n')
	repo = read_repo()
	if (func_str not in [x['func'] for x in repo]):
		with open(fp, "a") as ofile:
			ofile.write('#{}\n{}\n'.format(func_desc, func_str))

def _reload():
	init()
	get_ctx()['func_str'] = None
	func_str = ''
	if arg_has_key('-plot'):
		func_str = arg_get('-plot', 'sin(x)+cos(3*y)^2').decode('utf-8').encode('ascii', 'ignore'); func_str = ''.join(func_str.split());
		add_to_repo(func_str, arg_get('-desc', ''))
	else:
		repo = read_repo(def_repo_fp(), arg_get('-filter', '').split(','), arg_has('-with_scatter'))
		repol = sorted(repo, key=lambda el: el['i'])
		repoll = ['{} : {}'.format(el['desc'],el['func']) for el in repol]
		choices = printAndChoose(repoll, False, True, ' ')
		if len(choices):
			repo_choice = repol[choices[0]]
			func_str = repo_choice['func']
			get_ctx()['scatter'] = to_scatter(repo_choice['scatter'])
			get_ctx()['init_scatter'] = True
			get_ctx()['script'] = repo_choice['script']
			get_ctx()['init_script'] = True
	if arg_has_key('-cc'):
		get_ctx()['cc_str'] = arg_get('-cc', None)
	if len(func_str):
		post_init()
		is_console = sys.flags.interactive
		if (is_console):
			get_ctx()['style'] = 'console'
			plot(func_str if (get_ctx()['func_str'] is None) else get_ctx()['func_str'])
		else:
			get_ctx()['case_swap'] = arg_has('-case_swap')
			replt = True
			while(replt):
				replt = False
				plot(func_str if (get_ctx()['func_str'] is None) else get_ctx()['func_str'])
				replt = get_ctx()['replot']
				get_ctx()['replot'] = False

_reload()
