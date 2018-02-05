#!/bin/csh

foreach i ('actinium' 'aluminum' 'americium' 'antimony' 'argon' 'arsenic' 'astatine' 'barium' 'berkelium' 'beryllium' 'bismuth' 'boron' 'bromine' 'cadmium' 'calcium' 'californium' 'carbon' 'cerium' 'cesium' 'chlorine' 'chromium' 'cobalt' 'copper' 'curium' 'dysprosium' 'einsteinium' 'erbium' 'europium' 'fluorine' 'francium' 'gadolinium' 'gallium' 'germanium' 'gold' 'hafnium' 'helium' 'holmium' 'hydrogen' 'indium' 'iodine' 'iridium' 'iron' 'krypton' 'lanthanum' 'lead' 'lithium' 'lutetium' 'magnesium' 'manganese' 'mercury' 'molybdenum' 'neodymium' 'neon' 'neptunium' 'nickel' 'niobium' 'nitrogen' 'osmium' 'oxygen' 'palladium' 'phosphorus' 'platinum' 'plutonium' 'polonium' 'potassium' 'praseodymium' 'promethium' 'protactinium' 'radium' 'radon' 'rhenium' 'rhodium' 'rubidium' 'ruthenium' 'samarium' 'scandium' 'selenium' 'silicon' 'silver' 'sodium' 'strontium' 'sulfur' 'tantalum' 'technetium' 'tellurium' 'terbium' 'thallium' 'thorium' 'thulium' 'tin' 'titanium' 'tungsten' 'uranium' 'vanadium' 'xenon' 'ytterbium' 'yttrium' 'zinc' 'zirconium' )

    wget http://physics.nist.gov/PhysRefData/Handbook/Tables/findinglist.htm
    wget http://physics.nist.gov/PhysRefData/Handbook/element_name.htm
    wget http://physics.nist.gov/PhysRefData/Handbook/atomic_number.htm
    wget http://physics.nist.gov/PhysRefData/Handbook/periodictable.htm
    
    foreach j ( 1 2 3 4 5 6 7)
	wget "http://physics.nist.gov/PhysRefData/Handbook/Tables/"${i}"table"${j}".htm"
    end

end
