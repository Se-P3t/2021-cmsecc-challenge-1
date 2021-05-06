
# Result

环上序列的截位还原问题

## category 1

### level 1

### level 2

### level 3

### level 4

### level 5

### level 6

### level 7




## category 2

### level 1

```sh
 % time sage recover_coefficients__kernel.sage 2147483647 2 30 8 17 --category 2 --level 1 --verbose 1 --block-size 20 --check
SEED: 12928939533349146289

expect_vectors: 5

kernel rank: 2
kernel' rank: 2

find k: 7 0
checking z_i: maybe
checking c_i: True

================
CPU     102%
user    2.144
system  0.261
total   2.354
```

```json
{
    "modulus": 2147483647,
    "zbits": 17,
    "coefficients": [1596998372, 913674193],
    "initial_state": [25583676, 1935662022]
}
```

### level 2

```sh
 % time sage recover_coefficients__kernel.sage 2147483647 2 60 15 23 --category 2 --level 2 --verbose 1 --block-size 20 --check
SEED: 608460742630951138

expect_vectors: 5

kernel rank: 2
kernel' rank: 2

find k: -26204 -811
checking z_i: maybe
checking c_i: True

================
CPU     101%
user    2.714
system  0.237
total   2.902
```

```json
{
    "modulus": 2147483647,
    "zbits": 23,
    "coefficients": [423368878, 1375517413],
    "initial_state": [1968687461, 159779378]
}
```

### level 3

```sh
 % time sage recover_coefficients__kernel.sage 2147483647 3 68 17 21 --category 2 --level 3 --verbose 1 --block-size 20 --check
SEED: 15922487348339387766

expect_vectors: 5

kernel rank: 2
kernel' rank: 2

find k: -1610 -50
checking z_i: maybe
checking c_i: True

================
CPU     101%
user    2.914
system  0.296
total   3.160
```

```json
{
    "modulus": 2147483647,
    "zbits": 21,
    "coefficients": [233454232, 712694596, 1250324919],
    "initial_state": [1968687461, 159779378, 933973255]
}
```

### level 4

```sh
 % time sage recover_coefficients__kernel.sage 2147483647 4 95 25 21 --category 2 --level 4 --verbose 1 --block-size 30 --check 
SEED: 18135446848677888257

expect_vectors: 5

kernel rank: 2
kernel' rank: 2

find k: 1540 9
checking z_i: maybe
checking c_i: True

================
CPU     100%
user    7.603
system  0.266
total   7.828
```

```json
{
    "modulus": 2147483647,
    "zbits": 21,
    "coefficients": [1724886998, 287764120, 669496309, 29431304],
    "initial_state": [23214526, 63888791, 632526187, 515165230]
}
```

### level 5

```sh
 % time sage recover_coefficients__kernel.sage 2113941029 5 85 23 18 --category 2 --level 5 --verbose 1 --block-size 30 --check
SEED: 13512413672535046818

expect_vectors: 5

kernel rank: 2
kernel' rank: 2

find k: -24 0
checking z_i: maybe
checking c_i: True

================
CPU     100%
user    6.491
system  0.256
total   6.697
```

```json
{
    "modulus": 2113941029,
    "zbits": 18,
    "coefficients": [241592126, 1225700761, 270381722, 1809937814, 545364186],
    "initial_state": [1265879737, 1938499214, 100295411, 164486483, 782938]
}
```

### level 6

```sh
 % time sage recover_coefficients__kernel.sage 2140900439 8 90 20 11 --category 2 --level 6 --verbose 1 --block-size 20 --check
SEED: 2774189050861033585

expect_vectors: 6

kernel rank: 2

checking z_i: maybe
checking c_i: True

================
CPU     100%
user    4.257
system  0.153
total   4.403
```

```json
{
    "modulus": 2140900439,
    "zbits": 11,
    "coefficients": [1468898684, 429201201, 1438911747, 1343646518, 197478154, 1760674261, 1954064960, 1521596057],
    "initial_state": [2100894922, 1278709211, 1882769045, 822236456, 503871792, 639865048, 740172666, 1225091197]
}
```

### level 7

```sh
 % time sage recover_coefficients__kernel.sage 2086596509 10 110 26 11 --category 2 --level 7 --verbose 1 --block-size 20 --check
SEED: 4768417682194502211

expect_vectors: 6

kernel rank: 2

checking z_i: maybe
checking c_i: True

================
CPU     100%
user    9.136
system  0.159
total   9.289
```

```json
{
    "modulus": 2086596509,
    "zbits": 11,
    "coefficients": [1111873952, 1210156476, 822665602, 1221543376, 50425289, 1211476215, 1888560914, 1679919063, 1715756131, 141246785],
    "initial_state": [1800774673, 1719445571, 1277627869, 633482595, 1079260842, 2060264416, 573852671, 1245793554, 141816054, 891093089]
}
```

### level 8

```sh
 % time sage recover_coefficients__kernel.sage 2123058169 12 110 28 8 --category 2 --level 8 --verbose 1 --block-size 20 --check
SEED: 12270745752118493764

expect_vectors: 5

kernel rank: 2

checking z_i: maybe
checking c_i: True

================
CPU     100%
user    4.516
system  0.145
total   4.652
```

```json
{
    "modulus": 2123058169,
    "zbits": 8,
    "coefficients": [1380518532, 572802739, 397998604, 1517367287, 1838517468, 1894991963, 186761507, 1691926163, 917271042, 1840254211, 1520485994, 1544456793],
    "initial_state": [1407121795, 1964561578, 1522896740, 706799210, 1022864344, 1609032902, 945477963, 1469966024, 232347115, 2025154173, 653675533, 603475034]
}
```

### level 9

```sh
 % time sage recover_coefficients__kernel.sage 2147483647 14 128 32 8 --category 2 --level 9 --verbose 1 --block-size 32 --check
SEED: 8421578719646212792

expect_vectors: 5

kernel rank: 2

checking z_i: maybe
checking c_i: True

================
CPU     100%
user    9.820
system  0.190
total   10.004
```

```json
{
    "modulus": 2147483647,
    "zbits": 8,
    "coefficients": [755735009, 435105367, 1987422269, 141113323, 1831273687, 150474978, 1521781010, 88703098, 2128502238, 1314935750, 1897202874, 765777736, 1257457888, 851182418],
    "initial_state": [2043273873, 1844618476, 93566804, 974010131, 1813982016, 210537594, 1382311654, 397516522, 1569721206, 1648699957, 971748769, 1513210834, 1522349751, 1075882738]
}
```


## category 3

### level 1

```sh
 % time sage recover_modulus__kernel.sage 31 16 140 30 5 --category 3 --level 1 --verbose 1 --block-size 20 --check
SEED: 4630477721370439504

expect_vectors: 6

kernel rank: 2

checking z_i: maybe
finding modulus
 31 bits maybe
checking c_i: True

================
CPU     100%
user    1:05.23
system  0.334
total   1:05.57
```

```json
{
    "modulus": 2123847813,
    "zbits": 5,
    "coefficients": [2018484748, 1845189071, 1463275556, 1602150465, 194422107, 1025586312, 1724843525, 410393693, 106087782, 229852172, 1293380914, 235543842, 1642599451, 250756393, 239416449, 1593118903],
    "initial_state": [1422595032, 698794309, 2264898, 2123217555, 1515919515, 2048279701, 1227019859, 1625549320, 475257352, 1682624639, 1669847210, 510649846, 272012336, 1608958057, 1318317428, 1804116296]
}
```

### level 2

```sh
 % time sage recover_modulus__kernel.sage 31 16 190 40 10 --category 3 --level 2 --verbose 1 --block-size 30 --check
SEED: 5206466853635617849

expect_vectors: 6

kernel rank: 2

checking z_i: maybe
finding modulus
 31 bits maybe
checking c_i: True

================
CPU     99%
user    16:42.22
system  0.433
total   16:42.74
```

```json
{
    "modulus": 2146390813,
    "zbits": 10,
    "coefficients": [1709517653, 1473447434, 146866621, 1301261246, 1098198483, 1586144939, 880631859, 1190804449, 419206704, 377180855, 997067781, 668707083, 1991423106, 1200018419, 1124879071, 2081342702],
    "initial_state": [963105734, 789585151, 1238195227, 2028522939, 1124205863, 1618865668, 452174891, 77673612, 46901163, 1112184517, 2114462053, 259959215, 1976235589, 1517149832, 147549104, 73665604]
}
```
