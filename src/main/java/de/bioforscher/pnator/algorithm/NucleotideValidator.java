package de.bioforscher.pnator.algorithm;

import de.bioforscher.singa.structure.model.identifiers.LeafIdentifier;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * @author cl
 */
public class NucleotideValidator {

    private static final String[] namesToValidate = new String[]{"N1'", "C2'", "C3'", "N4'", "C7'", "C8'", "C5'", "C'"};

    private LeafIdentifier leafIdentifier;
    private HashMap<String, Boolean> replacedAtomNames;

    public NucleotideValidator(LeafIdentifier leafIdentifier) {
        this.leafIdentifier = leafIdentifier;
        this.replacedAtomNames = new HashMap<>();
        for (String name : namesToValidate) {
            replacedAtomNames.put(name, false);
        }
    }

    public void validate(String replacedName) {
        replacedAtomNames.put(replacedName, true);
    }

    public Set<String> getInvalidNames() {
        Set<String> invalidNames = new HashSet<>();
        for (Map.Entry<String, Boolean> entry : replacedAtomNames.entrySet()) {
            if (!entry.getValue()) {
                invalidNames.add(entry.getKey());
            }
        }
        return invalidNames;
    }

    public boolean isValid() {
        for (Map.Entry<String, Boolean> entry : replacedAtomNames.entrySet()) {
            if (!entry.getValue()) {
                return false;
            }
        }
        return true;
    }

}
