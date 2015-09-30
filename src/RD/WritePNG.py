#---------------------------------------------------------------------------
#
# WritePNG.py: writes compressed, true-color RGBA PNG files
#
# RGBA stands for "Red Green Blue Alpha", where alpha is the opacity level
#
# extracted from:
# http://stackoverflow.com/questions/902761/saving-a-numpy-array-as-an-image
#
# Original source code:
# https://developer.blender.org/diffusion/B/browse/master/release/bin/blender-thumbnailer.py$155
#

def write_png(buf, width, height):
    # by ideasman42, 2013-10-04, stackoverflow.com
    """ buf: must be bytes or a bytearray in py3, a regular string in py2. formatted RGBARGBA... """
    import zlib, struct

    # reverse the vertical line order and add null bytes at the start
    width_byte_4 = width * 4
    raw_data = b''.join(b'\x00' + buf[span:span + width_byte_4]
                        for span in range((height - 1) * width * 4, -1, - width_byte_4))

    def png_pack(png_tag, data):
        chunk_head = png_tag + data
        return (struct.pack("!I", len(data)) +
                chunk_head +
                struct.pack("!I", 0xFFFFFFFF & zlib.crc32(chunk_head)))

    return b''.join([
        b'\x89PNG\r\n\x1a\n',
        png_pack(b'IHDR', struct.pack("!2I5B", width, height, 8, 6, 0, 0, 0)),
        png_pack(b'IDAT', zlib.compress(raw_data, 9)),
        png_pack(b'IEND', b'')])

def test_write_png():
    # a red square:
    buf=b'\xFF\x00\x00\xFF'
    n=9
    imgsize=2**n # generate an image of size imgsize x imgsize pixels
    for i in range(2*n):
        buf = buf + buf
    print "len=", len(buf)/4

    # The data should be written directly to a file opened as binary, as in:
    data = write_png(buf, imgsize, imgsize)
    with open("my_image.png", 'wb') as fd:
        fd.write(data)

def saveAsPNG(array, filename):
    # by Evgeni Sergeev, 2014-01-10, stackoverflow.com
    import struct
    if any([len(row) != len(array[0]) for row in array]):
        raise ValueError, "Array should have elements of equal size"

                                #First row becomes top row of image.
    flat = []; map(flat.extend, reversed(array))
                                 #Big-endian, unsigned 32-byte integer.
    buf = b''.join([struct.pack('>I', ((0xffFFff & i32)<<8)|(i32>>24) )
                    for i32 in flat])   #Rotate from ARGB to RGBA.

    data = write_png(buf, len(array[0]), len(array))
    f = open(filename, 'wb')
    f.write(data)
    f.close()

def test_save_png():
    import numpy as np
    a = np.empty((2,2), np.uint32)
    a.fill(0xFF)
    r = np.empty((2,2), np.uint32)
    r[0,0] = 0xFF
    r[0,1] = 0xFF
    g = np.empty((2,2), np.uint32)
    g[0,1] = 0xFF
    b = np.empty((2,2), np.uint32)
    b[1,1] = 0xFF
    tot = (a << 24) | (r << 16) | (g << 8) | b
    print tot
    saveAsPNG(tot, 'test_grid.png')
    #saveAsPNG([[0xffff0000, 0xffFFFF00],
    #           [0xff00aa77, 0xff333333]], 'test_grid.png')

if __name__ == '__main__':
    test_save_png()

