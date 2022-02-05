// Print formatted content with a newline

#include "print.h"
#include "debug.h"
#include "noncopyable.h"
#include <fcntl.h>
#include <memory>
#include <unistd.h>
namespace mandelbrot {

using std::unique_ptr;

namespace {
struct Tee : public Noncopyable {
  FILE* f = 0;
  ~Tee() { if (f) fclose(f); }
};
}

static Tee& tee_fd() {
  static Tee t;
  return t;
}

void tee(const string& path) {
  auto& t = tee_fd();
  if (t.f) die("can't tee to two files");
  const int fd = open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0666);
  if (fd < 0) die("can't tee output to '%s': %s", path, strerror(errno));
  t.f = fdopen(fd, "w");
  if (!t.f) die("fdopen failed: %s", strerror(errno));
}

static inline void tee_print(const string& s) {
  if (FILE* f = tee_fd().f) {
    fprintf(f, "%s\n", s.c_str());
    fflush(f);
  }
}

void print() { print(""); }

void print(const string& s) {
  printf("%s\n", s.c_str());
  tee_print(s);
}

void print_error(const string& s) {
  fprintf(stderr, "%s\n", s.c_str());
  tee_print(s);
}

}  // namespace mandelbrot
