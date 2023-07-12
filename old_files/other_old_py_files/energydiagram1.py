#! /usr/bin/env python
#>import energydiagram as energyd


from energydiagram import ED
import matplotlib.pyplot as plt
diagram = ED()
diagram.add_level(0.0000000000,'0',)
diagram.add_level(0.1088966327,'1',)
diagram.add_level(0.1238192474,'2','last',)
diagram.add_level(0.1506515791,'3','last',)
diagram.add_level(0.1672213509,'4','last',)
diagram.add_level(0.1724090430,'5','last',)
diagram.add_level(0.1776359942,'6',)
diagram.add_level(0.1813476288,'7',)
diagram.add_level(0.1840561557,'8',)
diagram.add_level(0.1855497774,'9',)
diagram.add_level(0.1878920096,'10',)
diagram.add_level(0.1905931443,'11',)
diagram.add_level(0.1921082661,'12',)
diagram.add_level(0.2031759728,'13',)
diagram.add_level(0.2137368276,'14',)
diagram.add_level(0.2166105925,'15',)
diagram.add_level(0.2179305433,'16',)
diagram.add_level(0.2199858545,'17',)
diagram.add_level(0.2275083282,'18',)
diagram.add_level(0.2461181338,'19',)
diagram.add_level(0.2471716736,'20',)
diagram.add_level(0.2515222214,'21',)
diagram.add_level(0.2516920920,'22',)
diagram.add_level(0.2615077955,'23',)
diagram.add_level(0.2619076080,'24',)


diagram.plot(ylabel="Energy / $au$ ",show_IDs=True) # this is the default ylabel
plt.show()