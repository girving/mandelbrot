// Bottcher logo

#include "complex.h"
#include "expansion.h"
#include "series.h"
#include <png.h>
#include <array>
#include <random>

using namespace mandelbrot;
using std::array;
using std::mt19937;
using std::uniform_real_distribution;

typedef array<double,4> Color;
const Color clear = {0,0,0,0};

struct Box {
  Complex<double> lo, hi;

  Complex<double> center() const { return half(lo + hi); }
  Complex<double> shape() const { return hi - lo; }
};

// https://en.wikipedia.org/wiki/Alpha_compositing
static inline Color over(const Color a, const Color b) {
  Color r;
  const auto sa = a[3];
  const auto sb = (1 - sa)*b[3];
  const auto s = sa + sb;
  const auto inv_s = s ? 1 / s : 0;
  r[3] = s;
  for (int i = 0; i < 3; i++)
    r[i] = inv_s * (sa*a[i] + sb*b[i]);
  return r;
}

static inline Color average(span<const Color> cs) {
  Color sum = {0};
  for (const auto& c : cs) {
    sum[3] += c[3];
    for (int a = 0; a < 3; a++)
      sum[a] += c[a] * c[3];
  }
  Color mean;
  const auto scale = sum[3] ? 1 / sum[3] : 0;
  for (int a = 0; a < 3; a++)
    mean[a] = scale * sum[a];
  mean[3] = sum[3] / cs.size();
  return mean;
}

struct Canvas : public Noncopyable {
  Box box;
  double dx, inv_dx;
  int width, height, samples;
  Array<Complex<double>> points;
  Array<Color> colors;

  Canvas(const Box box_in, const int size, const int samples)
    : samples(samples) {
    // Adjust dimensions to have uniform aspect ratio
    const auto center = box_in.center();
    const auto shape = box_in.shape();
    dx = min(shape.r, shape.i) / size;
    inv_dx = 1 / dx;
    width = int(ceil(inv_dx * shape.r));
    height = int(ceil(inv_dx * shape.i));
    const auto h = half(dx) * Complex<double>(width, height);
    box = Box{center - h, center + h};
    print("canvas:");
    print("  samples %d", samples);
    print("  input: box %g %g, size %d", box_in.lo, box_in.hi, size);
    print("  final: box %g %g, size %d %d", box.lo, box.hi, width, height);

    // Compute samples for antialiasing
    Array<Complex<double>>(width * height * samples).swap(points);
    mt19937 mt;
    uniform_real_distribution<double> uniform;
    for (int i = 0; i < width; i++) {
      for (int j = 0; j < height; j++) {
        for (int s = 0; s < samples; s++) {
          const auto lo = box.lo + dx * Complex<double>(i, j);
          const auto x = lo.r + dx * uniform(mt);
          const auto y = lo.i + dx * uniform(mt);
          points[index(i,j,s)] = Complex<double>(x, y);
        }
      }
    }

    // Prepare for rendering!
    Array<Color>(width * height * samples).swap(colors);
  }

  int index(const int i, const int j, const int s) const {
    return (i*height + j)*samples + s;
  }

  // Draw an arbitrary shape
  template<class F> void render(const Box b, F&& f) const {
    const auto ilo = inv_dx * (b.lo - box.lo);
    const auto ihi = inv_dx * (b.hi - box.lo);
    const int i0 = max(0, int(floor(ilo.r)));
    const int j0 = max(0, int(floor(ilo.i)));
    const int i1 = min(width, 1+int(floor(ihi.r)));
    const int j1 = min(height, 1+int(floor(ihi.i)));
    for (int i = i0; i < i1; i++) {
      for (int j = j0; j < j1; j++) {
        for (int s = 0; s < samples; s++) {
          const int I = index(i,j,s);
          const auto z = points[I];
          auto& c = colors[I];
          c = over(f(z), c);
        }
      }
    }
  }

  // Draw a circle
  void circle(const Complex<double> center, const double radius, const Color color) const {
    const Complex<double> h(radius, radius);
    render(Box{center - h, center + h}, [center,radius,color](const Complex<double> z) {
      return sqr_abs(z - center) <= sqr(radius) ? color : clear;
    });
  }

  // Draw a Gaussian blob
  void gaussian(const Complex<double> center, const double radius, const Color color) const {
    const double scale = -0.5 / sqr(radius);
    const double bound = 5*radius;
    const Complex<double> h(bound, bound);
    render(Box{center - h, center + h}, [center,scale,color](const Complex<double> z) {
      const double s = exp(scale * sqr_abs(z - center));
      return array{color[0], color[1], color[2], s*color[3]};
    });
  }

  // Color the border
  void border(const Color color) const {
    slow_assert(width && height);
    for (int i = 0; i < width; i++)
      for (int s = 0; s < samples; s++)
        colors[index(i,0,s)] = colors[index(i,height-1,s)] = color;
    for (int j = 0; j < height; j++)
      for (int s = 0; s < samples; s++)
        colors[index(0,j,s)] = colors[index(width-1,j,s)] = color;
  }

  // Write to png
  void write(const string& path) const {
    FILE* f = fopen(path.c_str(), "wb");
    slow_assert(f, "failed to open %s for writing", path);
    auto png = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
    auto info = png_create_info_struct(png);
    png_init_io(png, f);
    png_set_IHDR(png, info, width, height, 8, PNG_COLOR_TYPE_RGB_ALPHA, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);
    const Array<uint8_t> row(width*4);
    for (int j = 0; j < height; j++) {
      for (int i = 0; i < width; i++) {
        const auto c = average(span<const Color>(&colors[index(i,j,0)], samples));
        for (int a = 0; a < 4; a++)
          row[i*4 + a] = max(0, min(255, int(rint(255 * c[a]))));
      }
      png_write_row(png, row.data());
    }
    png_destroy_write_struct(&png, &info);
    fclose(f);
  }
};

Complex<double> cis_tau(const double t) {
  const auto s = 2 * M_PI * t;
  return Complex<double>(cos(s), sin(s));
}

void logo() {
  typedef Expansion<2> E;
  typedef Complex<double> C;

  // Read f series
  const int max_k = 25;
  const auto [_, f_] = read_series<E>(format("exp2-11mar/f-k%d", max_k));
  const auto f = f_.view();

  // Render
  const int size = 256;
  const int samples = 256;
  const Canvas canvas(Box{{-2.01,-1.14},{.5,1.14}}, size, samples);
  //const Canvas canvas(Box{{.25,-.8},{2.02,.8}}, size, samples);
  const auto render = [f,&canvas](const int k, const double radius, const Color color) {
    // f[:2^k].astype(double)
    print("k %d", k);
    const int p = 1 << k;
    Series<double> fk(p);
    fk.set_counts(p, p);
    for (int i = 0; i < p; i++)
      fk[i] = double(f[i]);

    // Do an srfft to get point samples along the circle
    Array<C> fz(p/2);
    srfft<double>(fz, fk);

    // Render
    for (int i = 0; i < p/2; i++) {
      const auto z = cis_tau((i + 0.5) / p);
      const auto phi = z * fz[i];
      for (const auto c : {phi, conj(phi)})
        canvas.circle(c, radius, color);
    }
  };
  render(max_k, 0.001, Color{0,0,1,1});
  render(20, 0.002, Color{0,0,1,1});
  render(15, 0.003, Color{0,.7,.2,1});
  render(10, 0.006, Color{0,.2,1,1});
  render(7, 0.01, Color{0,1,.7,1});
  render(5, 0.02, Color{0,.9,.3,1});

  // Write to file
  canvas.write("logo.png");
}

int main(const int argc, const char** argv) {
  try {
    logo();
    return 0;
  } catch (const std::exception& e) {
    die(e.what());
  }
}
