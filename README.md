Some tools for encoding/decoding tiny vector images. See the 'shaders' folder
for shader functions.

<small>Note: this assumes you're running on a little-endian processor!</small>

To build:

    mkdir build && cd build
    cmake ..
    make

Using:

    # encoding an image with default settings, outputs a PNG at 'bar.png.vectex'
    ./vecimg -e foo.png bar.png.vectex
    # set block size, maps 16x16 blocks to each output pixel
    ./vecimg -e -b 16 foo.png bar.png.vectex
    # decoding an image (defaults to blocksize of 8)
    ./vecimg -d bar.png.vectex baz.png
    # decoding with explicit block size, mapping each input pixel to
    # 16x16 pixel output blocks
    ./vecimg -d -b 16 bar.png.vectex baz.png
