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
add_global_link_arguments('-lc++', language: 'cpp')

cpp = meson.get_compiler('cpp')
gmp = cpp.find_library('gmp', dirs: libpaths)
m = cpp.find_library('m')
flint = declare_dependency(link_args: ['-lflint'] + linkflags, dependencies: [gmp])
arb = declare_dependency(link_args: [darwin ? '-larb' : '-lflint-arb'] + linkflags, dependencies: [flint, m])
threads = dependency('threads')

cuda = dependency('cuda', version: '>=11.3', required: false)
if cuda.found()
  add_global_arguments('-I/usr/local/cuda/include', language: 'cpp')
  add_global_arguments('-x', 'cuda', '--cuda-gpu-arch=sm_80', language: 'cpp')
  add_global_arguments('-Wno-unknown-cuda-version', language: 'cpp')
endif

# Download tinyformat.h
tinyformat = custom_target('tinyformat', output: 'tinyformat.h',
                           command: ['curl', ('https://raw.githubusercontent.com/c42f/tinyformat/' +
                                     'aef402d85c1e8f9bf491b72570bfe8938ae26727/tinyformat.h'), '-o', '@OUTPUT@'])

known = declare_dependency(link_with: library('known', ['known.cc']))
base = declare_dependency(
    link_with: library('base', ['debug.cc', 'shutdown.cc', tinyformat]),
    sources: [tinyformat])
tests = declare_dependency(link_with: library('tests', ['tests.cc'], dependencies: [base, m]), dependencies: [base, m])
cutil = declare_dependency(
    link_with: library('cutil', ['cutil.cc'], dependencies: [base, cuda]),
    dependencies: [base, cuda])
arb_cc = declare_dependency(
    link_with: library('arbcc', ['acb_cc.cc', 'arb_cc.cc', 'arf_cc.cc', 'fmpq_cc.cc', 'poly.cc'],
                       dependencies: [base, arb]),
    dependencies: [base, arb])
arb_area = declare_dependency(
    link_with: library('arb_area', ['arb_area.cc'], dependencies: [arb_cc, known]),
    dependencies: [base, arb, arb_cc, known])

codelets = executable('codelets', ['codelets.cc'], dependencies: [base, m])
generated = [
    custom_target('gen_expansion', output: 'gen-expansion.h', command: [codelets, '@OUTPUT@']),
    custom_target('gen_mul_bases', output: 'gen-mul-bases.h', command: [codelets, '@OUTPUT@']),
]
area = declare_dependency(
    link_with: library('area', ['area.cc', 'expansion.cc', 'series.cc', 'fft.cc', 'nearest.cc'] + generated,
                       dependencies: [base, arb_cc, cutil, known, m]),
    sources: generated, dependencies: [arb_cc, cutil, known, m])

test('bit', executable('bit_test', ['bit_test.cc'], dependencies: [tests]))
test('expansion', executable('expansion_test', ['expansion_test.cc'], dependencies: [area, tests]))
test('arb_area', executable('arb_area_test', ['arb_area_test.cc'], dependencies: [arb_area, tests]))
test('poly', executable('poly_test', ['poly_test.cc'], dependencies: [arb_area, tests]))
test('fft', executable('fft_test', ['fft_test.cc'], dependencies: [area, tests]))
test('series', executable('series_test', ['series_test.cc'], dependencies: [area, tests]))
test('area', executable('area_test', ['area_test.cc'], dependencies: [area, arb_area, tests]))
if cuda.found()
  test('series_cuda', executable('series_cuda_test', ['series_cuda_test.cc'], dependencies: [area, cutil, tests]))
  test('cuda', executable('cuda_test', ['cuda_test.cc'], dependencies: [cutil, tests]))
  test('area_cuda', executable('area_cuda_test', ['area_cuda_test.cc'], dependencies: [area, arb_area, tests]))
endif

executable('mandelbrot', ['mandelbrot.cc'], dependencies: [area, arb_area])
