
# Result

环上序列的截位还原问题

## cpuinfo

```sh
 % cat /proc/cpuinfo | grep 'model name'
model name      : Intel(R) Core(TM) i7-8750H CPU @ 2.20GHz
```

## category 1

### level 1

```sh
 % time python recover_initial_state__embedding.py 2147483647 "257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768" 22 2 --category 1 --level 1 --verbose 1 --block-size 20      
SEED: 16949532060652715486

row 2/23:: True

solution: [257122259, 561754033, 340567466, 370127561, 1603890151, 789710361, 438276282, 1205614745, 929435387, 273101854, 1330188497, 1927005651, 70974738, 512638222, 655376420, 1799812840]

================
CPU     384%
user    2.044
system  1.548
total   0.935
```

### level 2

```sh
 % time python recover_initial_state__embedding.py 2147483647 "257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768" 23 3 --category 1 --level 2 --verbose 1 --block-size 20
SEED: 15864885015374467406

row 1/24:: True

solution: [720759561, 571519362, 1471239584, 2146374199, 1943406434, 1021387208, 911847040, 972284090, 1269846527, 2125700056, 582856260, 1326807684, 820744516, 83875554, 1033926054, 974403091]

================
CPU     375%
user    2.052
system  1.552
total   0.960
```

### level 3

```sh
 % time python recover_initial_state__embedding.py 2147483647 "257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768" 30 14 --category 1 --level 3 --verbose 1 --block-size 20
SEED: 10677130670124483019

row 1/31:: True

solution: [2056319062, 1175325906, 34344908, 300877541, 871395921, 953051611, 1276817066, 429446330, 716236050, 1498148665, 419137803, 983513800, 877501144, 162784581, 615817441, 1532450981]

================
CPU     356%
user    2.052
system  1.496
total   0.994
```

### level 4

```sh
 % time python recover_initial_state__embedding.py 2147483647 "257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768" 51 21 --category 1 --level 4 --verbose 1 --block-size 20
SEED: 6253750390576407687

row 5/52:: False
row 36/52:: False
row 39/52:: True

solution: [146123803, 1660954690, 553686861, 1592631770, 2039784960, 874444650, 1462760700, 1629573947, 927239148, 2020986341, 2134682761, 1440980008, 214415113, 823589071, 1178840115, 237668181]

================
CPU     327%
user    2.234
system  1.582
total   1.163
```

### level 5

```sh
 % time python recover_initial_state__embedding.py 2147483647 "257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768" 75 24 --category 1 --level 5 --verbose 1 --block-size 38
SEED: 8435949126609154252

row 1/76:: True

solution: [1683538635, 247852287, 1110454234, 2134379965, 524624671, 1135659173, 611030817, 917303344, 67431860, 844532212, 2121384016, 1630029172, 1311537205, 1289944654, 358213366, 2024311471]

================
CPU     206%
user    6.824
system  1.874
total   4.203
```

### level 6

```sh
 % time python recover_initial_state__embedding.py 2147483647 "257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768" 128 26 --category 1 --level 6 --verbose 2 --block-size 30 --threads 6 --sieve
SEED: 18336376286588014776

not found

Loaded file 'svpchallenge-129.txt'
gh = 657627253520767647744.000000, goal_r0/gh = 1.102500, r0/gh = 4.352077
  33: ↑ 33 ↓  5  T:   0.77867s, TT:   0.77870s, r0:2.26859e+21 r0/gh:   3.44967
  36: ↑ 36 ↓  1  T:   1.00746s, TT:   1.79837s, r0:2.26859e+21 r0/gh:   3.44967
  39: ↑ 39 ↓  4  T:   1.44775s, TT:   3.25843s, r0:2.26859e+21 r0/gh:   3.44967
  42: ↑ 42 ↓  1  T:   2.08314s, TT:   5.35389s, r0:2.03131e+21 r0/gh:   3.08884
  45: ↑ 45 ↓  6  T:   2.38020s, TT:   7.74711s, r0:2.03131e+21 r0/gh:   3.08884
  48: ↑ 48 ↓  8  T:   3.13874s, TT:  10.90063s, r0:1.98669e+21 r0/gh:   3.02100
  51: ↑ 51 ↓  1  T:   5.31238s, TT:  16.22622s, r0:1.98669e+21 r0/gh:   3.02100
  54: ↑ 54 ↓ 14  T:   5.03884s, TT:  21.27855s, r0:1.98669e+21 r0/gh:   3.02100
  57: ↑ 57 ↓ 22  T:   5.14259s, TT:  26.43470s, r0:1.85512e+21 r0/gh:   2.82092
  60: ↑ 60 ↓ 11  T:   6.19231s, TT:  32.63975s, r0:1.83267e+21 r0/gh:   2.78679
  63: ↑ 63 ↓ 23  T:   6.50527s, TT:  39.15742s, r0:1.62436e+21 r0/gh:   2.47003
  66: ↑ 66 ↓ 38  T:   7.18083s, TT:  46.35055s, r0:1.62436e+21 r0/gh:   2.47003
  69: ↑ 69 ↓ 39  T:  11.54635s, TT:  57.91050s, r0:1.34953e+21 r0/gh:   2.05212
  72: ↑ 72 ↓ 39  T:  18.55924s, TT:  76.48368s, r0:1.34953e+21 r0/gh:   2.05212
  75: ↑ 69       T:   4.03360s, TT:  80.53204s, r0:2.08450e+20 r0/gh:   0.31697
svp: norm 14437807996.7 ,hf 0.56300
'threads': 6,        :: n: 129, cputime 239.9787s, walltime: 80.5325s, flast: 54.00, |db|: 2^16.62
row 1/129:: True

solution: [512741665, 1096178369, 808807049, 608491186, 420891056, 1682771835, 358966452, 55989687, 890238631, 2137448551, 2058494244, 38743896, 96170410, 49379854, 1551639435, 1878181314]

================
CPU     288%
user    4:06.99
system  3.491
total   1:26.95
```

### level 7

for the multi-threaded version, the solution cannot be found after 10 threads run for about 12 hours

```sh
time python recover_initial_state__embedding.py 2147483647 "257,0,0,0,1048576,0,0,0,0,0,2097152,0,0,131072,0,32768" 150 27 --category 1 --level 7 --verbose 2 --block-size 30 --threads 10 --sieve
```

while sieving over GPU, a Tesla T4, we can get the results in 2 hours

POC at colab: https://colab.research.google.com/drive/1MwqKIxTzkqJMMD8vx3MhHY0l-8g976rm?usp=sharing

log: https://drive.google.com/file/d/1Y1B7usFfgONYfj78Hnt5POUDI5xMHEXT/view?usp=sharing

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
SEED: 10019631545210073147
USE_SUBS: False

find the kernel (rank 110)

expect_vectors: 6

find the kernel (rank 2)

kernel rank: 2

checking z_i: maybe
finding modulus
 31 bits maybe
checking c_i: True

================
CPU     278%
user    1:25.33
system  0.483
total   30.772
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
SEED: 2732523400034528389
USE_SUBS: False

find the kernel (rank 150)

expect_vectors: 6

find the kernel (rank 2)

kernel rank: 2

checking z_i: maybe
finding modulus
 31 bits maybe
checking c_i: True

================
CPU     318%
user    3:57.23
system  0.699
total   1:14.81
```

```json
{
    "modulus": 2146390813,
    "zbits": 10,
    "coefficients": [1709517653, 1473447434, 146866621, 1301261246, 1098198483, 1586144939, 880631859, 1190804449, 419206704, 377180855, 997067781, 668707083, 1991423106, 1200018419, 1124879071, 2081342702],
    "initial_state": [963105734, 789585151, 1238195227, 2028522939, 1124205863, 1618865668, 452174891, 77673612, 46901163, 1112184517, 2114462053, 259959215, 1976235589, 1517149832, 147549104, 73665604]
}
```

### level 3

```sh
 % time sage recover_modulus__kernel.sage 31 16 265 70 14 --category 3 --level 3 --verbose 2 --block-size 30 --threads 10 --sieve --timeout 300 --max-dim 80 --check
SEED: 11316827758390952469
USE_SUBS: False

find the kernel (rank 195)

Loaded file 'svpchallenge-195.txt'
gh = 165086.056394, goal_r0/gh = 0.000000, r0/gh = 0.909447
'threads': 10,       :: n: 195, cputime  0.0090s, walltime:  0.0090s, flast: -1.00, |db|: 2^0.00
expect_vectors: 5

find the kernel (rank 2)

kernel rank: 2

checking z_i: maybe
finding modulus
 31 bits maybe
checking c_i: True

================
CPU     815%
user    43:59.48
system  3.875
total   5:23.99
```

```json
{
    "modulus": 2147385873,
    "zbits": 14,
    "coefficients": [1433408125, 1187896517, 1493677871, 674581842, 873135996, 1326711093, 1726922573, 833159114, 897246300, 1821464147, 1671306427, 570111420, 268554940, 1034987293, 1263557877, 894951295],
    "initial_state": [1976346479, 1550017240, 380413913, 160326487, 379904755, 622685282, 377075770, 1527069890, 946605673, 2096385369, 1160042120, 619378670, 50681874, 1753705306, 1683069046, 1512092243]
}
```
