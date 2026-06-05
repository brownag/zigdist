# Makefile for zigdist package development

R ?= R
RSCRIPT ?= Rscript
R_LIBS ?= /home/andrew/.local/share/R/library
export R_LIBS

# Default target
all: document build check readme

.PHONY: help document test build check install clean readme vignettes fmt all

help:
	@echo "zigdist package development Makefile"
	@echo "Usage: make [target]"
	@echo "Targets:"
	@echo "  document   - Regenerate Rd documentation from roxygen comments"
	@echo "  test       - Run unit tests using tinytest"
	@echo "  build      - Build the source package tarball (R CMD build)"
	@echo "  check      - Run R CMD check on the latest built tarball"
	@echo "  install    - Install the package (R CMD INSTALL)"
	@echo "  clean      - Remove built tarballs, check logs, and temporary build files"
	@echo "  readme     - Render README.Rmd to README.md"
	@echo "  vignettes  - Build vignettes"
	@echo "  fmt        - Format Zig source code using 'zig fmt'"
	@echo "  all        - Run document, build, and check"

document:
	$(RSCRIPT) -e "roxygen2::roxygenise()"

test:
	$(RSCRIPT) -e "tinytest::run_test_dir('tests')"

build:
	$(R) CMD build .

check:
	@tarball=$$(ls -t zigdist_*.tar.gz 2>/dev/null | head -n 1); \
	if [ -z "$$tarball" ]; then \
		echo "No built tarball found. Running 'make build' first..."; \
		$(R) CMD build .; \
		tarball=$$(ls -t zigdist_*.tar.gz 2>/dev/null | head -n 1); \
	fi; \
	$(R) CMD check $$tarball

install:
	$(R) CMD INSTALL .

clean:
	rm -rf zigdist_*.tar.gz zigdist.Rcheck/

readme:
	$(RSCRIPT) -e "rmarkdown::render('README.Rmd')"

vignettes:
	$(RSCRIPT) -e "tools::buildVignettes(dir = '.')"

fmt:
	zig fmt src/zig/
