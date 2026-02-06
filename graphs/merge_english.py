import matplotlib.pyplot as plt
import numpy as np

# Set up larger font sizes and figure
plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 15,
    'ytick.labelsize': 15,
    'legend.fontsize': 20
})
categories = []
# english.50mb

#y_RLZsize = np.array([9505, 13523, 16694, 19860, 23386, 27583, 31835, 37466, 46611, 53609, 68472, 100000])
#y_reference = np.array([52248, 26536, 13012, 6142, 2808, 1230, 498, 169, 62, 20, 5, 0])

#xml.50mb

#y_RLZsize = np.array([4110, 6808, 9060, 11461, 14111, 22492, 36746, 54091, 82070, 100000])
#y_reference = np.array([46223, 20565, 8773, 3303, 940, 305, 77, 21, 4, 0])

# english.50mb in bits
#y_RLZsize = np.array([209110, 270460,317186, 357480, 397562, 413745, 413855, 449592, 466110, 428872, 410832, 300000])
#y_reference = np.array([417984, 212288, 104096, 49136, 22464, 9840, 3984, 1352, 496, 160, 40, 0])

#xml in bits
#y_RLZsize = np.array([469954000.0, 464932000.0, 454164000.0, 444667000.0, 446693000.0, 442434000.0, 445125000.0, 444800000.0, 446577000.0, 453981000.0, 454108000.0, 454131000.0, 467064000.0, 482862000.0, 480764000.0, 505906000.0, 487228000.0, 482568000.0, 568228000.0, 629146000.0])
#y_reference = np.array([151092000.0, 89374000.0, 52334200.0, 30082600.0, 16912000.0, 9284130.0, 4953640.0, 2573640.0, 1293490.0, 629296.0, 295920.0, 127112.0, 53352.0, 21000.0, 6928.0, 2264.0, 768.0, 200.0, 32.0, 0.0])

#english full

#Sarscov

#y_reference = np.array([1455540000.0, 717484000.0, 364359000.0, 181283000.0, 89038800.0, 43292700.0, 21437300.0, 10428300.0, 4958880.0, 2507910.0, 1306940.0, 743096.0, 438904.0, 266232.0, 181264.0, 125240.0, 90632.0, 66544.0, 48944.0, 35656.0, 24864.0, 17624.0, 12176.0, 8288.0, 5560.0, 3368.0, 2176.0, 1352.0, 784.0, 392.0, 160.0, 88.0, 32.0])
#y_RLZsize = np.array([1955770.0, 2671900.0, 3256150.0, 3794880.0, 4357180.0, 4944820.0, 5474880.0, 6421060.0, 8653620.0, 12284700.0, 112907000.0, 403713000.0, 706657000.0, 878380000.0, 1035170000.0, 1132900000.0, 1268470000.0, 1373660000.0, 1347330000.0, 1424220000.0, 1449790000.0, 1519250000.0, 1536400000.0, 1620080000.0, 1670340000.0, 1695550000.0, 1831840000.0, 1859550000.0, 1949920000.0, 1952360000.0, 2245990000.0, 2350050000.0, 2488200000.0])

y_reference = np.array([363884000.0, 179371000.0, 91089700.0, 45320800.0, 22259700.0, 10823200.0, 5359310.0, 2607080.0, 1239720.0, 626978.0, 326736.0, 185774.0, 109726.0, 66558.0, 45316.0, 31310.0, 22658.0, 16636.0, 12236.0, 8914.0, 6216.0, 4406.0, 3044.0, 2072.0, 1390.0, 842.0, 544.0, 338.0, 196.0, 98.0, 40.0, 22.0, 8.0])
y_RLZsize = np.array([1730110.0, 2351270.0, 2857440.0, 3320520.0, 3800950.0, 4285520.0, 4728300.0, 5525100.0, 7387240.0, 10442000.0, 93551700.0, 330311000.0, 569885000.0, 702704000.0, 820999000.0, 890134000.0, 996657000.0, 1079300000.0, 1036410000.0, 1095560000.0, 1101840000.0, 1154630000.0, 1152300000.0, 1215060000.0, 1234600000.0, 1233130000.0, 1332240000.0, 1328250000.0, 1364950000.0, 1301580000.0, 1453290000.0, 1468780000.0, 1421830000.0])
sums = y_RLZsize + y_reference
for s in range(0,y_reference.size):
    categories.append(s)



fig, ax = plt.subplots(figsize=(14, 8))

bars1 = ax.bar(categories, y_reference, label='Reference', 
               color='steelblue', edgecolor='black', linewidth=0.5)
bars2 = ax.bar(categories, y_RLZsize, bottom=y_reference, label='RLZ', 
               color='green', edgecolor='black', linewidth=0.5)

#threshold for english.50mb

#threshold_value = 18119  

#threshold for xml.50mb

#threshold_value = 9329
lz77 = 3457038
rlz_prefix = 21157760
ax.axhline(y=lz77, color='black', linewidth=3, linestyle='-', 
           label=f'z = {lz77:,}')

ax.axhline(y=rlz_prefix, color='red', linewidth=3, linestyle='-', 
           label=f'z̄ = {rlz_prefix:,}')

sumsChoice = float("inf")
ia_choice = 0

for ia, sumsPrev in enumerate(sums):
    if sumsPrev < sumsChoice:
        sumsChoice = sumsPrev
        ia_choice = ia

ax.annotate(
    str(sumsChoice),
    xy=(ia_choice, sumsChoice),
    xytext=(0, 5),          # 5 points upward
    textcoords="offset points",
    ha="center",
    fontsize=14,
    fontweight="bold"
)



#for i, total in enumerate(sums):
 #   ax.text(i, total + (max(sums)*0.01), f'{total:,}', 
 #           ha='center', va='bottom', fontsize=14, fontweight='bold')

#for i, (ref, rlz) in enumerate(zip(y_reference, y_RLZsize)):
    # Label for reference (bottom part)
 #   if ref > 0:
  #      ax.text(i, ref/2, f'{ref:,}', ha='center', va='center', 
   #             fontsize=14, color='white')
    
   # if rlz > 0:
    #    ax.text(i, ref + rlz/2, f'{rlz:,}', ha='center', va='center', 
     #           fontsize=14, color='white')


ax.set_title('SarsCoV2, first 400MB', fontsize=25, pad=20)
ax.set_xlabel('Iter.', fontsize=25, )
ax.set_ylabel('Bits', fontsize=25)

ax.legend(loc='lower left', framealpha=0.9)

ax.grid(True, axis='y', alpha=0.3, linestyle='--')

ax.set_ylim(0, max(sums) * 1.05)

ax.set_yscale('symlog', linthresh=1e6)
plt.tight_layout()
plt.show()