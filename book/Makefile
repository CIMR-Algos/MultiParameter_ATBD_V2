all:
	JUPYTER_BOOK_BUILD=true jupyter-book build .
pdf:
	JUPYTER_BOOK_BUILD=true jupyter-book build . --builder pdfhtml

pdfl:
	JUPYTER_BOOK_BUILD=true jupyter-book build . --builder pdflatex

clean:
	rm -rf ./_build
