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

struct Box {
  Complex<double> lo, hi;

  Complex<double> center() const { return half(lo + hi); }
  Complex<double> shape() const { return hi - lo; }
};

// https://en.wikipedia.org/wiki/Alpha_compositing
static inline Color over(const Color a, const Color b) {
  Color r;
  const auto sa = a[3];
  const auto sb = (1 - a[3])*b[3];
  r[3] = sa + sb;
  for (int i = 0; i < 3; i++)
    r[i] = (sa*a[i] + sb*b[i]) / r[3];
  return r;
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

  // Draw a Gaussian blob
  void render_gaussian(const Complex<double> center, const double radius, const Color color) const {
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
    const double scale = 255 / samples;
    const Array<uint8_t> row(width*4);
    for (int j = 0; j < height; j++) {
      for (int i = 0; i < width; i++) {
        Color sum{0};
        for (int s = 0; s < samples; s++) {
          const auto c = colors[(i*height + j)*samples + s];
          for (int a = 0; a < 4; a++)
            sum[a] += c[a];
        }
        for (int a = 0; a < 4; a++)
          row[i*4 + a] = max(0, min(255, int(rint(scale * sum[a]))));
      }
      png_write_row(png, row.data());
    }
    png_destroy_write_struct(&png, &info);
    fclose(f);
  }
};

void logo() {
  typedef Expansion<2> E;
  typedef Complex<double> C;

  // Read f series
  const int max_k = 20;
  const auto [_, f_] = read_series<E>(format("exp2-11mar/f-k%d", max_k));
  const auto f = f_.view();

  // Render
  const int size = 256;
  const int samples = 64;
  const Canvas canvas(Box{{.25,-.8},{2.02,.8}}, size, samples);
  const auto render_k = [f,&canvas](const int k, const double radius, const Color color) {
    // f[:2^k].astype(double)
    print("k %d", k);
    const int p = 1 << k;
    Series<double> fk(p);
    fk.set_counts(p, p);
    for (int i = 0; i < p; i++)
      fk[i] = double(f[i]);

    // Do an srfft to get point samples along the circle
    Array<C> zs(p);
    srfft<double>(zs, fk);

    // Render
    for (int i = 0; i < p; i++) {
      const auto phi = zs[i];
      for (const auto c : {phi, conj(phi)})
        canvas.render_gaussian(c, radius, color);
    }
  };
  render_k(20, 0.0015, Color{0,0,1,1});

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
