project('mandelbrot', 'cpp',
        version: '0.1',
        default_options: [
            'cpp_std=c++20',
            'warning_level=2',
            'werror=true',
        ])

# Library search paths
libpaths = []

# Are we on Mac?
darwin = build_machine.system() == 'darwin'
if darwin
  # Add Homebrew paths
  brew = run_command('brew', '--prefix', check: true).stdout().strip()
  add_global_arguments('-I' + brew + '/include', language: 'cpp')
  libpaths += [brew + '/lib']
endif

# Linker flags
linkflags = []
foreach p : libpaths
  linkflags += ['-L' + p]
endforeach

add_global_arguments('-stdlib=libc++', language: 'cpp')
add_global_link_arguments('-stdlib=libc++', language: 'cpp')  # Not the default on Linux
add_global_arguments('-Wsign-compare', language: 'cpp')
add_global_arguments('-Wno-gnu-statement-expression', language: 'cpp')
add_global_arguments('-Wno-unused-parameter', language: 'cpp')
add_global_arguments('-Wno-deprecated-copy', language: 'cpp')
add_global_arguments('-Wno-c++11-narrowing', language: 'cpp')
add_global_arguments('-Wno-unqualified-std-cast-call', language: 'cpp')

cpp = meson.get_compiler('cpp')
gmp = cpp.find_library('gmp', dirs: libpaths)
m = cpp.find_library('m')
flint = declare_dependency(link_args: ['-lflint'] + linkflags, dependencies: [gmp, m])
threads = dependency('threads')
openmp = darwin ? [] : [declare_dependency(link_args: linkflags, dependencies: [dependency('openmp')])]
crypto = dependency('libcrypto')
png = dependency('libpng')

cuda = dependency('cuda', version: '>=11.3', required: false)
if cuda.found()
  add_global_arguments('-I/usr/local/cuda/include', language: 'cpp')
  add_global_arguments('-x', 'cuda', '--cuda-gpu-arch=sm_80', '--cuda-gpu-arch=sm_90', language: 'cpp')
  add_global_arguments('-Wno-unknown-cuda-version', language: 'cpp')
  # CUDA ≥ 12 rejects libc++ on x86 unless asked; clang's CUDA support works with it
  add_global_arguments('-D_ALLOW_UNSUPPORTED_LIBCPP', language: 'cpp')
endif

# No implicit FMA contraction anywhere: it breaks expansion arithmetic, and without it CPU and GPU builds
# (and Mac and Linux) do identical arithmetic.  Clang defaults to -ffp-contract=fast with CUDA, and to on.
#   https://llvm.org/docs/CompileCudaWithLLVM.html#flags-that-control-numerical-code
add_global_arguments('-ffp-contract=off', language: 'cpp')

# Download tinyformat.h
tinyformat = custom_target('tinyformat', output: 'tinyformat.h',
                           command: ['curl', ('https://raw.githubusercontent.com/c42f/tinyformat/' +
                                     'aef402d85c1e8f9bf491b72570bfe8938ae26727/tinyformat.h'), '-o', '@OUTPUT@'])

# Download argparse.hpp
argparse = custom_target('argparse', output: 'argparse.hpp',
                         command: ['curl', ('https://raw.githubusercontent.com/p-ranav/argparse/' +
                                   '79df7808d2164b8d2186014c8de039c8b2dc00e7/include/argparse/argparse.hpp'),
                                   '-o', '@OUTPUT@'])

base_deps = openmp + (cuda.found() ? [cuda] : [])
base = declare_dependency(
    link_with: library('base', ['angles.cc', 'cutil.cc', 'debug.cc', 'loops.cc', 'numpy.cc', 'print.cc', 'shutdown.cc',
                                tinyformat, argparse], dependencies: base_deps),
    sources: [tinyformat, argparse], dependencies: base_deps)
tests = declare_dependency(
    link_with: library('tests', ['tests.cc'], dependencies: [base, m]),
    dependencies: [base, m])
known = declare_dependency(
    link_with: library('known', ['known.cc'], dependencies: [base]),
    dependencies: [base])
arb_cc = declare_dependency(
    link_with: library('arbcc', ['acb_cc.cc', 'arb_cc.cc', 'arf_cc.cc', 'fmpq_cc.cc', 'poly.cc'],
                       dependencies: [base, flint]),
    dependencies: [base, flint])
arb_area = declare_dependency(
    link_with: library('arb_area', ['arb_area.cc'], dependencies: [arb_cc, known]),
    dependencies: [base, flint, arb_cc, known])

codelets = executable('codelets', ['sig.cc', 'exp.cc', 'nearest.cc', 'codelets.cc'],
                      dependencies: [base, flint, arb_cc, m, crypto])
generated = [
    custom_target('gen_expansion', output: 'gen-expansion.h', command: [codelets, '@OUTPUT@']),
    custom_target('gen_mul_bases', output: 'gen-mul-bases.h', command: [codelets, '@OUTPUT@']),
    custom_target('gen_series_bases', output: 'gen-series-bases.h', command: [codelets, '@OUTPUT@']),
    custom_target('gen_butterflies', output: 'gen-butterflies.h', command: [codelets, '@OUTPUT@']),
]
area = declare_dependency(
    link_with: library('area', ['area.cc', 'expansion.cc', 'series.cc', 'fft.cc', 'nearest.cc', 'octaves.cc'] + generated,
                       dependencies: [base, arb_cc, known, m]),
    sources: generated, dependencies: [arb_cc, known, m])

test('angles', executable('angles_test', ['angles_test.cc'], dependencies: [tests]))
test('tests', executable('tests_test', ['tests_test.cc'], dependencies: [tests]))
test('escape', executable('escape_test', ['escape_test.cc', 'escape.cc'], dependencies: [tests]))
test('tree', executable('tree_test', ['tree_test.cc', 'tree.cc', 'engine.cc', 'escape.cc'], dependencies: [area, tests]))
test('bit', executable('bit_test', ['bit_test.cc'], dependencies: [tests]))
test('numpy', executable('numpy_test', ['numpy_test.cc'], dependencies: [tests]))
test('expansion', executable('expansion_test', ['expansion_test.cc'], dependencies: [area, tests]))
test('loops', executable('loops_test', ['loops_test.cc'], dependencies: [base, tests]))
test('arb_area', executable('arb_area_test', ['arb_area_test.cc'], dependencies: [arb_area, tests]))
test('lavaurs', executable('lavaurs_test', ['lavaurs_test.cc', 'lavaurs.cc'], dependencies: [area, tests]))
test('poly', executable('poly_test', ['poly_test.cc'], dependencies: [arb_area, tests]))
test('fft', executable('fft_test', ['fft_test.cc'], dependencies: [area, tests]))
test('octaves', executable('octaves_test', ['octaves_test.cc'], dependencies: [area, tests]))
test('series', executable('series_test', ['series_test.cc'], dependencies: [area, tests]))
test('area', executable('area_test', ['area_test.cc'], dependencies: [area, arb_area, tests]))
test('julia', executable('julia_test', ['julia_test.cc', 'julia.cc'], dependencies: [area, tests]))
test('bulb', executable('bulb_test', ['bulb_test.cc', 'bulb.cc', 'engine.cc'], dependencies: [area, tests]))
test('julia_tail', executable('julia_tail_test', ['julia_tail_test.cc', 'julia_tail.cc', 'engine.cc'], dependencies: [tests]))
if cuda.found()
  test('series_cuda', executable('series_cuda_test', ['series_cuda_test.cc'], dependencies: [area, tests]))
  test('cuda', executable('cuda_test', ['cuda_test.cc'], dependencies: [base, tests]))
  test('area_cuda', executable('area_cuda_test', ['area_cuda_test.cc'], dependencies: [area, arb_area, tests]))
endif

executable('mandelbrot', ['mandelbrot.cc'], dependencies: [area, arb_area])
executable('logo', ['logo.cc'], dependencies: [area, png])
executable('series-to-npy', ['series-to-npy.cc'], dependencies: [area])
executable('series-exp', ['series-exp.cc'], dependencies: [area])
executable('hyperbolic', ['hyperbolic.cc'], dependencies: [area])
executable('renorm', ['renorm.cc'], dependencies: [base])
executable('octscan', ['octscan.cc'], dependencies: [area])
executable('rootsum', ['rootsum.cc'], dependencies: [area])
executable('tuning', ['tuning.cc'], dependencies: [base])
executable('fattened', ['fattened.cc'], dependencies: [base])
executable('census', ['census.cc'], dependencies: [area])
executable('escape_area', ['escape_area.cc', 'escape.cc'], dependencies: [base])
executable('escape_tree', ['escape_tree.cc', 'tree.cc', 'engine.cc', 'escape.cc'], dependencies: [area])
executable('julia_area', ['julia_area.cc', 'julia.cc'], dependencies: [area])
executable('parabolic_band', ['analysis/parabolic_band.cc'])
executable('satellite_tail', ['analysis/satellite_tail.cc'])
executable('parabolic_tail', ['analysis/parabolic_tail.cc', 'julia_tail.cc', 'engine.cc'], dependencies: [base])
executable('root_escape', ['analysis/root_escape.cc'])
executable('multiplier_profile', ['analysis/multiplier_profile.cc'])
executable('bulb_batch', ['analysis/bulb_batch.cc', 'bulb.cc', 'engine.cc'], dependencies: [area])
executable('bulb_areas', ['analysis/bulb_areas.cc'], dependencies: [area])
executable('satellite_tree', ['analysis/satellite_tree.cc'])
executable('family_owner', ['analysis/family_owner.cc', 'engine.cc'], dependencies: [base])
executable('limb_families', ['analysis/limb_families.cc'], dependencies: [base, area])
executable('lavaurs_area', ['analysis/lavaurs/lavaurs_area.cc', 'lavaurs.cc'], dependencies: [area])
executable('glavaurs', ['analysis/lavaurs/glavaurs_main.cc', 'glavaurs.cc', 'glavaurs_gpu.cc'], dependencies: [area])
executable('orbit_bench', ['orbit_bench.cc', 'engine.cc'], dependencies: [base])
executable('orbit_census', ['orbit_census.cc', 'engine.cc', 'escape.cc'], dependencies: [base])
executable('orbit_regimes', ['orbit_regimes.cc', 'engine.cc', 'escape.cc'], dependencies: [base])
