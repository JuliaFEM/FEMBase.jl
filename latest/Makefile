all: $(wildcard *.md)

%.md: %.ipynb
	jupyter nbconvert --to=markdown $<
