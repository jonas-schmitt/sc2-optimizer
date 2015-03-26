- TODO Combat.txt & Units.txt: Protoss: Oracle, Tempest, MothershipCore
- TODO Combat.txt & Units.txt: Zerg: SwarmHost, Viper
- TODO Combat.txt & Units.txt: Terran: Hellbat, WidowMine

Hier sind ein paar Erläuterungen zu dem Aufbau der Daten mit etwaigen TODOs:

Units:
	-Die ersten 5 Spalten werd ich wohl nicht erklären müssen, deswegen ignorier ich die jetzt.
	-In der Spalte BuildFrom stehen die Buildings/Units die die entsprechenden Units bauen können. Da man Buildings upgraden kann, diese aber trotzdem weiterhin die vorherigen Units bauen können, reicht es, wenn eines der genannten Dinge vorhanden ist.
	-Die Spalte Requirements stellen die zusätzlichen Requirements dar, zu den Buildings/Units, mit denen die Units produziert werden. Diese bleiben allerdings nach dem Herstellungsprozess erhalten.
	-Weiterhin gibt es die Vanish_Req Spalte. Dort werden die Units aufgelistet, die benötigt werden, bei der Erzeugung der neuen Unit allerdings verschwinden. (z.B. bei Protoss die Archon, und natürlich bei den Zerg häufig anzutreffen)


Buildings:
	-Hier sind die Dateien genauso aufgelistet, jedoch ist die Reihenfolge zwischen der Supply und BuildTime Spalte umgekehrt. Beim einlesen also aufpassen!
	-Weiterhin gibt es keine BuildFrom-Spalte.
	-Wenn bei den Requirements die Buildings durch ein '/' getrennt sind, bedeutet dass, dass nur eines der durch '/' getrennten Buildings notwendig ist. Dadurch werden auch die Upgrades von bestimmten Buildings abgedeckt.


Combat:
	- Hier befinden sich die Units mit ihren Kampfwerten
	- TODO: Bei den Zerg der Baneling. Ich hab keine DPS gefunden, nur einen Attackwert. Da der Baneling aber kein Cooldown hat, hab ich einfach diesen Attackwert als DPS genommen. Vielleicht einfach nochmal drüber schauen
	- TODO: Bei den Protoss der Archon. Der unterschied bei den Archons sind die Kosten, je nachdem aus welchen Templarn die erzeugt wurden. Ich habe einfach die mittleren Kosten genommen. Die anderen Statuswerte sind bei den Archons sonst gleich, deswegen hab ich da keine weitere Unterscheidung gemacht, aber evtl noch diskussionsbedarf

