package de.bioforscher.pnator.algorithm;

import de.bioforscher.singa.chemistry.descriptive.elements.Element;
import de.bioforscher.singa.chemistry.parser.pdb.structures.StructureParser;
import de.bioforscher.singa.chemistry.parser.pdb.structures.StructureParserOptions;
import de.bioforscher.singa.chemistry.parser.pdb.structures.StructureWriter;
import de.bioforscher.singa.chemistry.physical.atoms.Atom;
import de.bioforscher.singa.chemistry.physical.atoms.AtomName;
import de.bioforscher.singa.chemistry.physical.atoms.RegularAtom;
import de.bioforscher.singa.chemistry.physical.leaves.Nucleotide;
import de.bioforscher.singa.chemistry.physical.model.Structure;
import de.bioforscher.singa.mathematics.vectors.Vector3D;
import de.bioforscher.singa.mathematics.vectors.Vectors3D;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import static de.bioforscher.pnator.algorithm.AtomNameEquivalents.*;
import static de.bioforscher.singa.chemistry.descriptive.elements.ElementProvider.*;

/**
 * @author Alexander Eisold
 * @version 1.0
 */

public class PNAGenerator {

    private static final Logger logger = LoggerFactory.getLogger(PNAGenerator.class);

    private static int addedAtomIndex = 1000;
    private static final double C_O_DOUBLE_BOND_DISTANCE = 1.21;

    private static int backboneFailCount = 0;

    public static void main(String[] args) {

        /*
         * TODO: decide whether structure is dna or rna or hybrid
         *
         */

/*
        Structure structure = StructureParser.online()
                .pdbIdentifier("1BNA").everything().setOptions(StructureParserOptions.withSettings(StructureParserOptions.Setting.OMIT_HYDROGENS))
                .parse();
*/
         Structure structure = StructureParser.local().inputStream(Thread.currentThread().getContextClassLoader()
         .getResourceAsStream("structure_examples/model0001.pdb")).allModels().everything().setOptions(StructureParserOptions.withSettings(StructureParserOptions.Setting.OMIT_HYDROGENS)).parse();
         logger.info("Parsing structure {}.", structure.getPdbIdentifier());
        convertToPNAStructure(structure);

        try {
            StructureWriter.writeBranchSubstructure(structure.getFirstModel(), Paths.get("/tmp/test5.pdb"));
        } catch (IOException e) {
            e.printStackTrace();
        }

        /*
         StructureViewer.colorScheme = ColorScheme.BY_ELEMENT;
         StructureViewer.structure = structure;
         Application.launch(StructureViewer.class);
         */

    }

    public static Structure convertToPNAStructure(Structure structure) {

        structure.getAllChains().forEach(chain -> {

            List<Nucleotide> nucleotides = chain.getNucleotides();
            logger.info("Collected {} nucleotides for chain {}.", nucleotides.size(), chain.getIdentifier());

            if (nucleotides.isEmpty()) {
                logger.info("Chain {} is no nucleosid, skipping.", chain.getIdentifier());
            } else {
                nucleotides.forEach((Nucleotide nucleotide) -> {

                    Atom oxygenTow = PNAGenerator.calculateMissingAtoms(nucleotide.getAtomByName(AtomName.getAtomNameFromString("C1'")), nucleotide.getAtomByName(AtomName.getAtomNameFromString("C3'")), nucleotide.getAtomByName(AtomName.getAtomNameFromString("C2'")), false, "O7'", OXYGEN);

                    nucleotide.addNode(oxygenTow);
                    nucleotide.addEdgeBetween(nucleotide.getAtomByName(AtomName.getAtomNameFromString("C3'")), oxygenTow);

                    Optional<Atom> firstPhosphateOptional = FIRST_BACKBONE_PHOSPHATE.getAtomFrom(nucleotide);
                    Optional<Atom> secondPhosphateOptional = SECOND_BACKBONE_PHOSPHATE.getAtomFrom(nucleotide);
                    Optional<Atom> backbonePosphateOptional = BACKBONE_PHOSPHATE.getAtomFrom(nucleotide);

                    //RNA specific atoms
                    Optional<Atom> backboneOxygenTwoPrimeOptional = OXYGEN_TWO_PRIME.getAtomFrom(nucleotide);
                    Optional<Atom> backboneHydrogenOxygenTwoPrimeOptional = HYDROGEN_OXYGEN_TWO_PRIME.getAtomFrom(nucleotide);




                    if (firstPhosphateOptional.isPresent() && secondPhosphateOptional.isPresent() &&
                            backbonePosphateOptional.isPresent()) {

                        Atom oxygenOne = PNAGenerator.calculateMissingAtoms(firstPhosphateOptional.get(),
                                secondPhosphateOptional.get(), backbonePosphateOptional.get(), true,
                                "O1'", OXYGEN);

                        nucleotide.addNode(oxygenOne);
                        nucleotide.addEdgeBetween(backbonePosphateOptional.get(), oxygenOne);

                    } else {
                        backboneFailCount++;
                        if (backboneFailCount == 1) {
                            logger.warn("Could not calculate backbone for nucleotide {}.", nucleotide);
                        } else {
                            throw new InvalidInputStructure("Missing atoms to calculate backbone.");
                        }
                    }

                    // remove obsolete atoms
                    firstPhosphateOptional.ifPresent(nucleotide::removeNode);
                    secondPhosphateOptional.ifPresent(nucleotide::removeNode);


                    // remove obsolete RNA specific atoms
                    backboneOxygenTwoPrimeOptional.ifPresent(nucleotide::removeNode);
                    backboneHydrogenOxygenTwoPrimeOptional.ifPresent(nucleotide::removeNode);


                    chain.removeNode(nucleotide.getAtomByName(AtomName.O4Pr));

                    nucleotide.getAllAtoms().forEach(PNAGenerator::convertAtom);

                });
            }

            backboneFailCount = 0;
        });
        return structure;

    }


    private static void convertAtom(Atom an) {

        String atomName = an.getAtomNameString();



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
