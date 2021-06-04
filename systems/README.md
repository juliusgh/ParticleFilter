# Dynamische Systeme
Dieser Ordner enthält die Matlab-Implementierung von verschiedenen dynamischen Systemen, an welchen die Regelungskonzepte simulativ erprobt wurden.

Folgende Systeme sind verfügbar:
- `System.m`: Die Basisklasse für alle dynamischen Systeme, enthält ein 2-dimensionales System zum Testen der Partikelfilter
- `IntegratorSystem.m`: Ein simpler Doppelintegrator
- `PendulumSystem.m`: System eines Doppelpendels, wobei der Zustand direkt gemessen wird
- `PendulumSystem2.m`: System eines Doppelpendels, wobei die Position und Geschwindigkeit des äußeren Pendelstabs gemessen wird
- `PendulumSystem3.m`: System eines Doppelpendels, wobei die Positionen der beiden Pendelstäbe gemessen werden
- `ManipulatorSystem.m`: System eines ebenen Manipulators, wobei die Positionen der beiden Stabenden gemessen werden
- `ManipulatorSystemDirect.m`: System eines ebenen Manipulators, wobei der Zustand direkt gemessen wird

Julius Herb (st160887@stud.uni-stuttgart.de)