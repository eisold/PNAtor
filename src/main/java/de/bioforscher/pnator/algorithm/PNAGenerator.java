package de.bioforscher.pnator.algorithm;

import de.bioforscher.singa.chemistry.descriptive.elements.Element;
import de.bioforscher.singa.chemistry.descriptive.elements.ElementProvider;
import de.bioforscher.singa.chemistry.parser.pdb.structures.StructureParser;
import de.bioforscher.singa.chemistry.parser.pdb.structures.StructureWriter;
import de.bioforscher.singa.chemistry.physical.atoms.Atom;
import de.bioforscher.singa.chemistry.physical.atoms.AtomName;
import de.bioforscher.singa.chemistry.physical.atoms.RegularAtom;
import de.bioforscher.singa.chemistry.physical.leaves.Nucleotide;
import de.bioforscher.singa.chemistry.physical.model.StructuralEntityFilter;
import de.bioforscher.singa.chemistry.physical.model.Structure;
import de.bioforscher.singa.javafx.viewer.ColorScheme;
import de.bioforscher.singa.javafx.viewer.StructureViewer;
import de.bioforscher.singa.mathematics.vectors.Vector3D;
import de.bioforscher.singa.mathematics.vectors.Vectors3D;
import javafx.application.Application;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.lang.reflect.Modifier;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static de.bioforscher.singa.chemistry.descriptive.elements.ElementProvider.*;

/**
 * @author Alexander Eisold
 * @version 1.0
 */

public class PNAGenerator {

    private static final Logger logger = LoggerFactory.getLogger(PNAGenerator.class);

    private static int addedAtomIndex = 1000;
    private static final double C_O_DOUBLE_BOND_DISTANCE = 1.21;

    public static void main(String[] args) {

        Structure structure = StructureParser.online()
                .pdbIdentifier("1BNA")
                .parse();

        logger.info("Parsing structure {}.", structure.getPdbIdentifier());

        List<Nucleotide> nucleotides = structure.getAllLeaves().stream()
                .filter(StructuralEntityFilter.LeafFilter.isNucleotide())
                .map(leaf -> (Nucleotide) leaf)
                .collect(Collectors.toList());

        logger.info("Collected {} nucleotides.", nucleotides.size());

        nucleotides
                .forEach((Nucleotide nucleotide) -> {

                    Atom oxygenTow = PNAGenerator.calculateMissingAtoms(nucleotide.getAtomByName(AtomName.getAtomNameFromString("C1'")), nucleotide.getAtomByName(AtomName.getAtomNameFromString("C3'")), nucleotide.getAtomByName(AtomName.getAtomNameFromString("C2'")), false, "O7'", OXYGEN);

                    nucleotide.addNode(oxygenTow);
                    nucleotide.addEdgeBetween(nucleotide.getAtomByName(AtomName.getAtomNameFromString("C3'")), oxygenTow);
                    structure.getFirstModel().get().removeNode(nucleotide.getAtomByName(AtomName.O4Pr));


                    Atom atomOP1 = null;
                    if (nucleotide.containsAtomWithName(AtomName.OP1)) {
                        atomOP1 = nucleotide.getAtomByName(AtomName.OP1);
                    }
                    Atom atomOP2 = null;
                    if (nucleotide.containsAtomWithName(AtomName.OP2)) {
                        atomOP2 = nucleotide.getAtomByName(AtomName.OP2);
                    }
                    Atom atomP = null;
                    if (nucleotide.containsAtomWithName(AtomName.P)) {
                        atomP = nucleotide.getAtomByName(AtomName.P);
                    }
                    if (atomOP1 != null && atomOP2 != null && atomP != null) {


                        Atom oxygenOne = PNAGenerator.calculateMissingAtoms(atomOP1, atomOP2, atomP, true, "O1'", OXYGEN);

                        nucleotide.addNode(oxygenOne);
                        nucleotide.addEdgeBetween(atomP, oxygenOne);

                        structure.getFirstModel().get().removeNode(atomOP1);
                        structure.getFirstModel().get().removeNode(atomOP2);


                    } else {
                        logger.warn("Could not calculate backbone for nucleotide {}.", nucleotide);
                    }


                    nucleotide.getAllAtoms().forEach(PNAGenerator::convertAtom);

                });


        try {
            StructureWriter.writeBranchSubstructure(structure.getFirstModel().get(), Paths.get("/tmp/test.pdb"));
        } catch (IOException e) {
            e.printStackTrace();
        }
/**
 StructureViewer.colorScheme = ColorScheme.BY_ELEMENT;
 StructureViewer.structure = structure;
 Application.launch(StructureViewer.class);
 **/

    }


    private static void convertAtom(Atom an) {

        switch (an.getAtomNameString()) {

            case "O5'":
                an.setElement(NITROGEN);
                an.setAtomNameString("N1'");
                logger.trace("Replacing Atom {} through {}.", "O5'", an.getAtomNameString());
                break;
            case "C5'":

                an.setElement(CARBON);
                an.setAtomNameString("C2'");
                logger.trace("Replacing Atom {} through {}.", "C5'", an.getAtomNameString());

                break;

            case "C4'":

                an.setElement(CARBON);
                an.setAtomNameString("C3'");
                logger.trace("Replacing Atom {} through {}.", "C4'", an.getAtomNameString());

                break;

            case "C3'":

                an.setElement(NITROGEN);
                an.setAtomNameString("N4'");
                logger.trace("Replacing Atom {} through {}.", "C3'", an.getAtomNameString());

                break;

            case "C2'":

                an.setElement(CARBON);
                an.setAtomNameString("C7'");
                logger.trace("Replacing Atom {} through {}.", "C2'", an.getAtomNameString());

                break;

            case "C1'":

                an.setElement(CARBON);
                an.setAtomNameString("C8'");
                logger.trace("Replacing Atom {} through {}.", "C1'", an.getAtomNameString());

                break;

            case "O3'":

                an.setElement(CARBON);
                an.setAtomNameString("C5'");
                logger.trace("Replacing Atom {} through {}.", "O3'", an.getAtomNameString());

                break;

            case "P":

                an.setElement(CARBON);
                an.setAtomNameString("C'");
                logger.trace("Replacing Atom {} through {}.", "P", an.getAtomNameString());

                break;

            default:
                break;
        }
    }

    /**
     * @pram a first Atom for the calculation of the centroid between Atom a and Atom b
     * @pram b second Atom for the calculation of the centroid between Atom a and Atom b
     * @pram c an Atom that is bounded with the calculated missing Atom
     * @pram forward an boolean that defined the direction of the new normalized vector (true is positive and false negate the vector)
     * @pram name an String of the missing atom
     * @pram element an Element equals the missing atom
     */

    public static Atom calculateMissingAtoms(Atom a, Atom b, Atom c, Boolean forward, String name, Element element) {

        Vector3D positionA = a.getPosition();
        Vector3D positionB = b.getPosition();

        Vector3D centroid = Vectors3D.getCentroid(Arrays.asList(positionA, positionB));
        logger.trace("Calculated centroid between {} and {} as {}.", a, b, centroid);

        Vector3D missingAtomPosition = null;
        if (forward) {
            missingAtomPosition = c.getPosition().add(centroid.subtract(c.getPosition()).normalize().multiply(C_O_DOUBLE_BOND_DISTANCE));

        } else {
            missingAtomPosition = c.getPosition().add(centroid.subtract(c.getPosition()).normalize().additivelyInvert().multiply(C_O_DOUBLE_BOND_DISTANCE));

        }

        Atom missingAtom = new RegularAtom(addedAtomIndex++, element, name, missingAtomPosition);

        return missingAtom;
    }
}
