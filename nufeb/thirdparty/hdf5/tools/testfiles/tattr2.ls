Opened "tattr2.h5" with sec2 driver.
dset                     Dataset {2/2}
    Attribute: array {2}
        Type:      [3] 32-bit little-endian integer
        Data:  [1,2,3], [4,5,6]
    Attribute: array2D {3, 2}
        Type:      [3] 32-bit little-endian integer
        Data:
            (0,0) [1,2,3], [4,5,6],
            (1,0) [7,8,9], [10,11,12],
            (2,0) [13,14,15], [16,17,18]
    Attribute: array3D {4, 3, 2}
        Type:      [3] 32-bit little-endian integer
        Data:
            (0,0,0) [1,2,3], [4,5,6],
            (0,1,0) [7,8,9], [10,11,12],
            (0,2,0) [13,14,15], [16,17,18],
            (1,0,0) [19,20,21], [22,23,24],
            (1,1,0) [25,26,27], [28,29,30],
            (1,2,0) [31,32,33], [34,35,36],
            (2,0,0) [37,38,39], [40,41,42],
            (2,1,0) [43,44,45], [46,47,48],
            (2,2,0) [49,50,51], [52,53,54],
            (3,0,0) [55,56,57], [58,59,60],
            (3,1,0) [61,62,63], [64,65,66],
            (3,2,0) [67,68,69], [70,71,72]
    Attribute: bitfield {2}
        Type:      8-bit bitfield
        Data:  0x01, 0x02
    Attribute: bitfield2D {3, 2}
        Type:      8-bit bitfield
        Data:
            (0,0) 0x01, 0x02,
            (1,0) 0x03, 0x04,
            (2,0) 0x05, 0x06
    Attribute: bitfield3D {4, 3, 2}
        Type:      8-bit bitfield
        Data:
            (0,0,0) 0x01, 0x02,
            (0,1,0) 0x03, 0x04,
            (0,2,0) 0x05, 0x06,
            (1,0,0) 0x07, 0x08,
            (1,1,0) 0x09, 0x0a,
            (1,2,0) 0x0b, 0x0c,
            (2,0,0) 0x0d, 0x0e,
            (2,1,0) 0x0f, 0x10,
            (2,2,0) 0x11, 0x12,
            (3,0,0) 0x13, 0x14,
            (3,1,0) 0x15, 0x16,
            (3,2,0) 0x17, 0x18
    Attribute: compound {2}
        Type:      struct {
                   "a"                +0    8-bit integer
                   "b"                +4    IEEE 64-bit little-endian float
               } 12 bytes
        Data:  {1, 2}, {3, 4}
    Attribute: compound2D {3, 2}
        Type:      struct {
                   "a"                +0    8-bit integer
                   "b"                +4    IEEE 64-bit little-endian float
               } 12 bytes
        Data:
            (0,0) {1, 2}, {3, 4},
            (1,0) {5, 6}, {7, 8},
            (2,0) {9, 10}, {11, 12}
    Attribute: compound3D {4, 3, 2}
        Type:      struct {
                   "a"                +0    8-bit integer
                   "b"                +4    IEEE 64-bit little-endian float
               } 12 bytes
        Data:
            (0,0,0) {1, 2}, {3, 4},
            (0,1,0) {5, 6}, {7, 8},
            (0,2,0) {9, 10}, {11, 12},
            (1,0,0) {13, 14}, {15, 16},
            (1,1,0) {17, 18}, {19, 20},
            (1,2,0) {21, 22}, {23, 24},
            (2,0,0) {25, 26}, {27, 28},
            (2,1,0) {29, 30}, {31, 32},
            (2,2,0) {33, 34}, {35, 36},
            (3,0,0) {37, 38}, {39, 40},
            (3,1,0) {41, 42}, {43, 44},
            (3,2,0) {45, 46}, {47, 48}
    Attribute: enum {2}
        Type:      enum 32-bit little-endian integer {
                   RED              = 0
                   GREEN            = 1
               }
        Data:  RED, RED
    Attribute: enum2D {3, 2}
        Type:      enum 32-bit little-endian integer {
                   RED              = 0
                   GREEN            = 1
               }
        Data:
            (0,0) RED, RED,
            (1,0) RED, RED,
            (2,0) RED, RED
    Attribute: enum3D {4, 3, 2}
        Type:      enum 32-bit little-endian integer {
                   RED              = 0
                   GREEN            = 1
               }
        Data:
            (0,0,0) RED, RED,
            (0,1,0) RED, RED,
            (0,2,0) RED, RED,
            (1,0,0) RED, RED,
            (1,1,0) RED, RED,
            (1,2,0) RED, RED,
            (2,0,0) RED, RED,
            (2,1,0) RED, RED,
            (2,2,0) RED, RED,
            (3,0,0) RED, RED,
            (3,1,0) RED, RED,
            (3,2,0) RED, RED
    Attribute: float {2}
        Type:      IEEE 32-bit little-endian float
        Data:  1, 2
    Attribute: float2D {3, 2}
        Type:      IEEE 32-bit little-endian float
        Data:
            (0,0) 1, 2,
            (1,0) 3, 4,
            (2,0) 5, 6
    Attribute: float3D {4, 3, 2}
        Type:      IEEE 32-bit little-endian float
        Data:
            (0,0,0) 1, 2,
            (0,1,0) 3, 4,
            (0,2,0) 5, 6,
            (1,0,0) 7, 8,
            (1,1,0) 9, 10,
            (1,2,0) 11, 12,
            (2,0,0) 13, 14,
            (2,1,0) 15, 16,
            (2,2,0) 17, 18,
            (3,0,0) 19, 20,
            (3,1,0) 21, 22,
            (3,2,0) 23, 24
    Attribute: integer {2}
        Type:      32-bit little-endian integer
        Data:  1, 2
    Attribute: integer2D {3, 2}
        Type:      32-bit little-endian integer
        Data:
            (0,0) 1, 2,
            (1,0) 3, 4,
            (2,0) 5, 6
    Attribute: integer3D {4, 3, 2}
        Type:      32-bit little-endian integer
        Data:
            (0,0,0) 1, 2,
            (0,1,0) 3, 4,
            (0,2,0) 5, 6,
            (1,0,0) 7, 8,
            (1,1,0) 9, 10,
            (1,2,0) 11, 12,
            (2,0,0) 13, 14,
            (2,1,0) 15, 16,
            (2,2,0) 17, 18,
            (3,0,0) 19, 20,
            (3,1,0) 21, 22,
            (3,2,0) 23, 24
    Attribute: opaque {2}
        Type:      1-byte opaque type
               (tag = "1-byte opaque type")
        Data:  0x01, 0x02
    Attribute: opaque2D {3, 2}
        Type:      1-byte opaque type
               (tag = "1-byte opaque type")
        Data:
            (0,0) 0x01, 0x02,
            (1,0) 0x03, 0x04,
            (2,0) 0x05, 0x06
    Attribute: opaque3D {4, 3, 2}
        Type:      1-byte opaque type
               (tag = "1-byte opaque type")
        Data:
            (0,0,0) 0x01, 0x02,
            (0,1,0) 0x03, 0x04,
            (0,2,0) 0x05, 0x06,
            (1,0,0) 0x07, 0x08,
            (1,1,0) 0x09, 0x0a,
            (1,2,0) 0x0b, 0x0c,
            (2,0,0) 0x0d, 0x0e,
            (2,1,0) 0x0f, 0x10,
            (2,2,0) 0x11, 0x12,
            (3,0,0) 0x13, 0x14,
            (3,1,0) 0x15, 0x16,
            (3,2,0) 0x17, 0x18
    Attribute: reference {2}
        Type:      object reference
        Data:  DATASET-1:976, DATASET-1:976
    Attribute: reference2D {3, 2}
        Type:      object reference
        Data:
            (0,0) DATASET-1:976, DATASET-1:976,
            (1,0) DATASET-1:976, DATASET-1:976,
            (2,0) DATASET-1:976, DATASET-1:976
    Attribute: reference3D {4, 3, 2}
        Type:      object reference
        Data:
            (0,0,0) DATASET-1:976, DATASET-1:976,
            (0,1,0) DATASET-1:976, DATASET-1:976,
            (0,2,0) DATASET-1:976, DATASET-1:976,
            (1,0,0) DATASET-1:976, DATASET-1:976,
            (1,1,0) DATASET-1:976, DATASET-1:976,
            (1,2,0) DATASET-1:976, DATASET-1:976,
            (2,0,0) DATASET-1:976, DATASET-1:976,
            (2,1,0) DATASET-1:976, DATASET-1:976,
            (2,2,0) DATASET-1:976, DATASET-1:976,
            (3,0,0) DATASET-1:976, DATASET-1:976,
            (3,1,0) DATASET-1:976, DATASET-1:976,
            (3,2,0) DATASET-1:976, DATASET-1:976
    Attribute: string {2}
        Type:      2-byte null-terminated ASCII string
        Data:  "ab", "de"
    Attribute: string2D {3, 2}
        Type:      2-byte null-terminated ASCII string
        Data:
            (0,0) "ab", "cd",
            (1,0) "ef", "gh",
            (2,0) "ij", "kl"
    Attribute: string3D {4, 3, 2}
        Type:      2-byte null-terminated ASCII string
        Data:
            (0,0,0) "ab", "cd",
            (0,1,0) "ef", "gh",
            (0,2,0) "ij", "kl",
            (1,0,0) "mn", "pq",
            (1,1,0) "rs", "tu",
            (1,2,0) "vw", "xz",
            (2,0,0) "AB", "CD",
            (2,1,0) "EF", "GH",
            (2,2,0) "IJ", "KL",
            (3,0,0) "MN", "PQ",
            (3,1,0) "RS", "TU",
            (3,2,0) "VW", "XZ"
    Attribute: vlen {2}
        Type:      variable length of
                   32-bit little-endian integer
        Data:  (1), (2,3)
    Attribute: vlen2D {3, 2}
        Type:      variable length of
                   32-bit little-endian integer
        Data:
            (0,0) (0), (1),
            (1,0) (2,3), (4,5),
            (2,0) (6,7,8), (9,10,11)
    Attribute: vlen3D {4, 3, 2}
        Type:      variable length of
                   32-bit little-endian integer
        Data:
            (0,0,0) (0), (1),
            (0,1,0) (2), (3),
            (0,2,0) (4), (5),
            (1,0,0) (6,7), (8,9),
            (1,1,0) (10,11), (12,13),
            (1,2,0) (14,15), (16,17),
            (2,0,0) (18,19,20), (21,22,23),
            (2,1,0) (24,25,26), (27,28,29),
            (2,2,0) (30,31,32), (33,34,35),
            (3,0,0) (36,37,38,39), (40,41,42,43),
            (3,1,0) (44,45,46,47), (48,49,50,51),
            (3,2,0) (52,53,54,55), (56,57,58,59)
    Location:  1:976
    Links:     1
    Storage:   8 logical bytes, 0 allocated bytes
    Type:      32-bit little-endian integer
g1                       Group
    Attribute: array {2}
        Type:      [3] 32-bit little-endian integer
        Data:  [1,2,3], [4,5,6]
    Attribute: array2D {3, 2}
        Type:      [3] 32-bit little-endian integer
        Data:
            (0,0) [1,2,3], [4,5,6],
            (1,0) [7,8,9], [10,11,12],
            (2,0) [13,14,15], [16,17,18]
    Attribute: array3D {4, 3, 2}
        Type:      [3] 32-bit little-endian integer
        Data:
            (0,0,0) [1,2,3], [4,5,6],
            (0,1,0) [7,8,9], [10,11,12],
            (0,2,0) [13,14,15], [16,17,18],
            (1,0,0) [19,20,21], [22,23,24],
            (1,1,0) [25,26,27], [28,29,30],
            (1,2,0) [31,32,33], [34,35,36],
            (2,0,0) [37,38,39], [40,41,42],
            (2,1,0) [43,44,45], [46,47,48],
            (2,2,0) [49,50,51], [52,53,54],
            (3,0,0) [55,56,57], [58,59,60],
            (3,1,0) [61,62,63], [64,65,66],
            (3,2,0) [67,68,69], [70,71,72]
    Attribute: bitfield {2}
        Type:      8-bit bitfield
        Data:  0x01, 0x02
    Attribute: bitfield2D {3, 2}
        Type:      8-bit bitfield
        Data:
            (0,0) 0x01, 0x02,
            (1,0) 0x03, 0x04,
            (2,0) 0x05, 0x06
    Attribute: bitfield3D {4, 3, 2}
        Type:      8-bit bitfield
        Data:
            (0,0,0) 0x01, 0x02,
            (0,1,0) 0x03, 0x04,
            (0,2,0) 0x05, 0x06,
            (1,0,0) 0x07, 0x08,
            (1,1,0) 0x09, 0x0a,
            (1,2,0) 0x0b, 0x0c,
            (2,0,0) 0x0d, 0x0e,
            (2,1,0) 0x0f, 0x10,
            (2,2,0) 0x11, 0x12,
            (3,0,0) 0x13, 0x14,
            (3,1,0) 0x15, 0x16,
            (3,2,0) 0x17, 0x18
    Attribute: compound {2}
        Type:      struct {
                   "a"                +0    8-bit integer
                   "b"                +4    IEEE 64-bit little-endian float
               } 12 bytes
        Data:  {1, 2}, {3, 4}
    Attribute: compound2D {3, 2}
        Type:      struct {
                   "a"                +0    8-bit integer
                   "b"                +4    IEEE 64-bit little-endian float
               } 12 bytes
        Data:
            (0,0) {1, 2}, {3, 4},
            (1,0) {5, 6}, {7, 8},
            (2,0) {9, 10}, {11, 12}
    Attribute: compound3D {4, 3, 2}
        Type:      struct {
                   "a"                +0    8-bit integer
                   "b"                +4    IEEE 64-bit little-endian float
               } 12 bytes
        Data:
            (0,0,0) {1, 2}, {3, 4},
            (0,1,0) {5, 6}, {7, 8},
            (0,2,0) {9, 10}, {11, 12},
            (1,0,0) {13, 14}, {15, 16},
            (1,1,0) {17, 18}, {19, 20},
            (1,2,0) {21, 22}, {23, 24},
            (2,0,0) {25, 26}, {27, 28},
            (2,1,0) {29, 30}, {31, 32},
            (2,2,0) {33, 34}, {35, 36},
            (3,0,0) {37, 38}, {39, 40},
            (3,1,0) {41, 42}, {43, 44},
            (3,2,0) {45, 46}, {47, 48}
    Attribute: enum {2}
        Type:      enum 32-bit little-endian integer {
                   RED              = 0
                   GREEN            = 1
               }
        Data:  RED, RED
    Attribute: enum2D {3, 2}
        Type:      enum 32-bit little-endian integer {
                   RED              = 0
                   GREEN            = 1
               }
        Data:
            (0,0) RED, RED,
            (1,0) RED, RED,
            (2,0) RED, RED
    Attribute: enum3D {4, 3, 2}
        Type:      enum 32-bit little-endian integer {
                   RED              = 0
                   GREEN            = 1
               }
        Data:
            (0,0,0) RED, RED,
            (0,1,0) RED, RED,
            (0,2,0) RED, RED,
            (1,0,0) RED, RED,
            (1,1,0) RED, RED,
            (1,2,0) RED, RED,
            (2,0,0) RED, RED,
            (2,1,0) RED, RED,
            (2,2,0) RED, RED,
            (3,0,0) RED, RED,
            (3,1,0) RED, RED,
            (3,2,0) RED, RED
    Attribute: float {2}
        Type:      IEEE 32-bit little-endian float
        Data:  1, 2
    Attribute: float2D {3, 2}
        Type:      IEEE 32-bit little-endian float
        Data:
            (0,0) 1, 2,
            (1,0) 3, 4,
            (2,0) 5, 6
    Attribute: float3D {4, 3, 2}
        Type:      IEEE 32-bit little-endian float
        Data:
            (0,0,0) 1, 2,
            (0,1,0) 3, 4,
            (0,2,0) 5, 6,
            (1,0,0) 7, 8,
            (1,1,0) 9, 10,
            (1,2,0) 11, 12,
            (2,0,0) 13, 14,
            (2,1,0) 15, 16,
            (2,2,0) 17, 18,
            (3,0,0) 19, 20,
            (3,1,0) 21, 22,
            (3,2,0) 23, 24
    Attribute: integer {2}
        Type:      32-bit little-endian integer
        Data:  1, 2
    Attribute: integer2D {3, 2}
        Type:      32-bit little-endian integer
        Data:
            (0,0) 1, 2,
            (1,0) 3, 4,
            (2,0) 5, 6
    Attribute: integer3D {4, 3, 2}
        Type:      32-bit little-endian integer
        Data:
            (0,0,0) 1, 2,
            (0,1,0) 3, 4,
            (0,2,0) 5, 6,
            (1,0,0) 7, 8,
            (1,1,0) 9, 10,
            (1,2,0) 11, 12,
            (2,0,0) 13, 14,
            (2,1,0) 15, 16,
            (2,2,0) 17, 18,
            (3,0,0) 19, 20,
            (3,1,0) 21, 22,
            (3,2,0) 23, 24
    Attribute: opaque {2}
        Type:      1-byte opaque type
               (tag = "1-byte opaque type")
        Data:  0x01, 0x02
    Attribute: opaque2D {3, 2}
        Type:      1-byte opaque type
               (tag = "1-byte opaque type")
        Data:
            (0,0) 0x01, 0x02,
            (1,0) 0x03, 0x04,
            (2,0) 0x05, 0x06
    Attribute: opaque3D {4, 3, 2}
        Type:      1-byte opaque type
               (tag = "1-byte opaque type")
        Data:
            (0,0,0) 0x01, 0x02,
            (0,1,0) 0x03, 0x04,
            (0,2,0) 0x05, 0x06,
            (1,0,0) 0x07, 0x08,
            (1,1,0) 0x09, 0x0a,
            (1,2,0) 0x0b, 0x0c,
            (2,0,0) 0x0d, 0x0e,
            (2,1,0) 0x0f, 0x10,
            (2,2,0) 0x11, 0x12,
            (3,0,0) 0x13, 0x14,
            (3,1,0) 0x15, 0x16,
            (3,2,0) 0x17, 0x18
    Attribute: string {2}
        Type:      2-byte null-terminated ASCII string
        Data:  "ab", "de"
    Attribute: string2D {3, 2}
        Type:      2-byte null-terminated ASCII string
        Data:
            (0,0) "ab", "cd",
            (1,0) "ef", "gh",
            (2,0) "ij", "kl"
    Attribute: string3D {4, 3, 2}
        Type:      2-byte null-terminated ASCII string
        Data:
            (0,0,0) "ab", "cd",
            (0,1,0) "ef", "gh",
            (0,2,0) "ij", "kl",
            (1,0,0) "mn", "pq",
            (1,1,0) "rs", "tu",
            (1,2,0) "vw", "xz",
            (2,0,0) "AB", "CD",
            (2,1,0) "EF", "GH",
            (2,2,0) "IJ", "KL",
            (3,0,0) "MN", "PQ",
            (3,1,0) "RS", "TU",
            (3,2,0) "VW", "XZ"
    Attribute: vlen {2}
        Type:      variable length of
                   32-bit little-endian integer
        Data:  (1), (2,3)
    Attribute: vlen2D {3, 2}
        Type:      variable length of
                   32-bit little-endian integer
        Data:
            (0,0) (0), (1),
            (1,0) (2,3), (4,5),
            (2,0) (6,7,8), (9,10,11)
    Attribute: vlen3D {4, 3, 2}
        Type:      variable length of
                   32-bit little-endian integer
        Data:
            (0,0,0) (0), (1),
            (0,1,0) (2), (3),
            (0,2,0) (4), (5),
            (1,0,0) (6,7), (8,9),
            (1,1,0) (10,11), (12,13),
            (1,2,0) (14,15), (16,17),
            (2,0,0) (18,19,20), (21,22,23),
            (2,1,0) (24,25,26), (27,28,29),
            (2,2,0) (30,31,32), (33,34,35),
            (3,0,0) (36,37,38,39), (40,41,42,43),
            (3,1,0) (44,45,46,47), (48,49,50,51),
            (3,2,0) (52,53,54,55), (56,57,58,59)
    Location:  1:2176
    Links:     1
g2                       Group
    Location:  1:2824
    Links:     1
