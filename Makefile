all: area

area: area.c
	$(CC) -o $@ -g -Wall -Werror -O2 -I$(HOME)/homebrew/include $^ -L$(HOME)/homebrew/lib -larb

.PHONY: clean
clean:
	rm -f area
