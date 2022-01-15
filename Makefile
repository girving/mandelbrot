all: area

area: area.c known.h
	$(CC) -o $@ -g -Wall -Werror -O2 -I$(HOME)/homebrew/include $< -L$(HOME)/homebrew/lib -larb -lflint

.PHONY: clean
clean:
	rm -f area
