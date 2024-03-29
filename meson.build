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
  add_global_arguments('-x', 'cuda', '--cuda-gpu-arch=sm_80', language: 'cpp')
  add_global_arguments('-Wno-unknown-cuda-version', language: 'cpp')

  # Clang defaults to -ffp-contract=fast with CUDA, which breaks expansion arithmetic
  #   https://llvm.org/docs/CompileCudaWithLLVM.html#flags-that-control-numerical-code
  add_global_arguments('-ffp-contract=off', language: 'cpp')
endif

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
    link_with: library('base', ['cutil.cc', 'debug.cc', 'loops.cc', 'numpy.cc', 'print.cc', 'shutdown.cc',
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
    link_with: library('area', ['area.cc', 'expansion.cc', 'series.cc', 'fft.cc', 'nearest.cc'] + generated,
                       dependencies: [base, arb_cc, known, m]),
    sources: generated, dependencies: [arb_cc, known, m])

test('bit', executable('bit_test', ['bit_test.cc'], dependencies: [tests]))
test('expansion', executable('expansion_test', ['expansion_test.cc'], dependencies: [area, tests]))
test('loops', executable('loops_test', ['loops_test.cc'], dependencies: [base, tests]))
test('arb_area', executable('arb_area_test', ['arb_area_test.cc'], dependencies: [arb_area, tests]))
test('poly', executable('poly_test', ['poly_test.cc'], dependencies: [arb_area, tests]))
test('fft', executable('fft_test', ['fft_test.cc'], dependencies: [area, tests]))
test('series', executable('series_test', ['series_test.cc'], dependencies: [area, tests]))
test('area', executable('area_test', ['area_test.cc'], dependencies: [area, arb_area, tests]))
if cuda.found()
  test('series_cuda', executable('series_cuda_test', ['series_cuda_test.cc'], dependencies: [area, tests]))
  test('cuda', executable('cuda_test', ['cuda_test.cc'], dependencies: [base, tests]))
  test('area_cuda', executable('area_cuda_test', ['area_cuda_test.cc'], dependencies: [area, arb_area, tests]))
endif

executable('mandelbrot', ['mandelbrot.cc'], dependencies: [area, arb_area])
executable('logo', ['logo.cc'], dependencies: [area, png])
executable('series-to-npy', ['series-to-npy.cc'], dependencies: [area])
executable('series-exp', ['series-exp.cc'], dependencies: [area])
