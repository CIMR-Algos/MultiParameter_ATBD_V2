html:
	JUPYTER_BOOK_BUILD=true jupyter-book build . -v

pdf:
	JUPYTER_BOOK_BUILD=true jupyter-book build . --builder pdfhtml -v 

pdfl:
	JUPYTER_BOOK_BUILD=true jupyter-book build . --builder pdflatex -v
	pdftk CIMR_DEVALGO_WP1_MultiParameterRetrievalATBD_D0_v2_frontpage.pdf _build/latex/atbd.pdf cat output CIMR_DEVALGO_WP1_MultiParameterRetrievalATBD_D0_v2.pdf

purge:
	rm -rf ./_build ./.ipynb_checkpoints ./figures

clean:
	jupyter-book clean .
