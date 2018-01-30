package de.bioforscher.pnator.algorithm;

import de.bioforscher.singa.structure.model.interfaces.Atom;
import de.bioforscher.singa.structure.model.interfaces.LeafSubstructure;
import de.bioforscher.singa.structure.model.oak.AtomName;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Optional;
import java.util.Set;

/**
 * @author cl
 */
public enum AtomNameEquivalents {

    FIRST_BACKBONE_PHOSPHATE("OP1", "O1P"),
    SECOND_BACKBONE_PHOSPHATE("OP2", "O2P"),
    BACKBONE_PHOSPHATE("P"),
    BACKBONE_OXYGEN_FIVE_PRIME("O5'"),
    BACKBONE_OXYGEN_THREE_PRIME("O3'"),
    BACKBONE_OXYGEN_FOUR_PRIME("O4'"),
    BACKBONE_CARBON_FIVE_PRIME("C5'"),
    BACKBONE_CARBON_FOUR_PRIME("C4'"),
    BACKBONE_CARBON_THREE_PRIME("C3'"),
    BACKBONE_CARBON_TWO_PRIME("C2'"),
    BACKBONE_CARBON_ONE_PRIME("C1'"),
    //RNA specific atoms
    BACKBONE_OXYGEN_TWO_PRIME("O2'");

    private Set<String> equivalentAtomNames;

    AtomNameEquivalents(String ... atomNames) {
        this.equivalentAtomNames = new HashSet<>();
        this.equivalentAtomNames.addAll(Arrays.asList(atomNames));
    }

    public Set<String> getEquivalentAtomNames() {
        return equivalentAtomNames;
    }

    public Optional<Atom> getAtomFrom(LeafSubstructure<?> leafSubstructure) {
        for (String equivalentAtomName : this.equivalentAtomNames) {
            AtomName atomName = AtomName.getAtomNameFromString(equivalentAtomName);
            return leafSubstructure.getAtomByName(atomName.getName());
        }
        return Optional.empty();
    }

    

}
