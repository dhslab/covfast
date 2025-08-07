# Makefile for the covfast program
# This Makefile will automatically download and compile the required version of htslib.

# --- Configuration ---
# Compiler and flags
CC = gcc
CFLAGS = -Wall -O3 -march=native
# Add the local htslib directory to the include path
CPPFLAGS = -I./htslib
# Link against the static htslib library directly, and its dependencies
LDFLAGS =
LDLIBS = ./htslib/libhts.a -lpthread -lz -lm -lbz2 -llzma -lcurl -lcrypto

# Program name
TARGET = covfast

# htslib configuration
HTSLIB_VERSION = 1.21
HTSLIB_URL = https://github.com/samtools/htslib/releases/download/$(HTSLIB_VERSION)/htslib-$(HTSLIB_VERSION).tar.bz2
HTSLIB_DIR = htslib-$(HTSLIB_VERSION)
HTSLIB_ARCHIVE = $(HTSLIB_DIR).tar.bz2
HTSLIB_LIB = ./htslib/libhts.a

# --- Targets ---

# Default target: build the covfast executable
all: $(TARGET)

# Rule to build the covfast executable
# This depends on the htslib library being built first.
# By explicitly providing the .a file, we ensure static linking.
$(TARGET): covfast.c $(HTSLIB_LIB)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(LDFLAGS) $(LDLIBS)

# Rule to build htslib
# This creates a static library 'libhts.a' in the ./htslib directory
$(HTSLIB_LIB):
	@echo "--- Downloading and building htslib $(HTSLIB_VERSION) ---"
	# Download the archive if it doesn't exist
	@if [ ! -f "$(HTSLIB_ARCHIVE)" ]; then \
		wget -O $(HTSLIB_ARCHIVE) $(HTSLIB_URL); \
	fi
	# Extract the archive
	tar -xjf $(HTSLIB_ARCHIVE)
	# Configure and build htslib
	cd $(HTSLIB_DIR) && ./configure --prefix=$(realpath .) && make
	# Create a symlink for a consistent path
	@ln -sf $(HTSLIB_DIR) htslib
	@echo "--- htslib build complete ---"

# Rule to install the executable to /usr/local/bin
install: $(TARGET)
	@echo "Installing $(TARGET) to /usr/local/bin..."
	install -m 0755 $(TARGET) /usr/local/bin/

# Rule to clean up build artifacts
clean:
	@echo "Cleaning up..."
	rm -f $(TARGET) *.o
	rm -rf $(HTSLIB_DIR) htslib $(HTSLIB_ARCHIVE)

.PHONY: all clean install
