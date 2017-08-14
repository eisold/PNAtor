package de.bioforscher.pnator.algorithm;

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


public class PNAGenerator {

    private static final Logger logger = LoggerFactory.getLogger(PNAGenerator.class);

    private static int addedAtomIndex = 1000;

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
                        Vector3D positionOP1 = atomOP1.getPosition();
                        Vector3D positionOP2 = atomOP2.getPosition();

                        Vector3D centroid = Vectors3D.getCentroid(Arrays.asList(positionOP1, positionOP2));
                        logger.trace("Calculated centroid between {} and {} as {}.",positionOP1, positionOP2,
                                centroid);


                        Vector3D newPosition = atomP.getPosition().add(centroid.subtract(atomP.getPosition()).normalize().multiply(1.0));
                        Atom atomOone = new RegularAtom(addedAtomIndex++, ElementProvider.OXYGEN,
                                "O1'", newPosition);
                        nucleotide.addNode(atomOone);
                        nucleotide.addEdgeBetween(atomP,atomOone);

                        System.out.println(positionOP1.angleTo(atomOone.getPosition()));

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

        StructureViewer.colorScheme = ColorScheme.BY_ELEMENT;
        StructureViewer.structure = structure;
        Application.launch(StructureViewer.class);


    }



    private static void convertAtom(Atom an) {

        switch (an.getAtomNameString()){

            case "O5'":
                an.setElement(ElementProvider.NITROGEN);
                an.setAtomNameString("N1'");
                logger.trace("Replacing Atom {} through {}.", "O5'",an.getAtomNameString() );
                break;
            case "C5'":

                an.setElement(ElementProvider.CARBON);
                an.setAtomNameString("C2'");
                logger.trace("Replacing Atom {} through {}.", "C5'",an.getAtomNameString() );

                break;

            case "C4'":

                an.setElement(ElementProvider.CARBON);
                an.setAtomNameString("C3'");
                logger.trace("Replacing Atom {} through {}.", "C4'",an.getAtomNameString() );

                break;

            case "C3'":

                an.setElement(ElementProvider.NITROGEN);
                an.setAtomNameString("N4'");
                logger.trace("Replacing Atom {} through {}.", "C3'",an.getAtomNameString() );

                break;

            case "C2'":

                an.setElement(ElementProvider.CARBON);
                an.setAtomNameString("C7'");
                logger.trace("Replacing Atom {} through {}.", "C2'",an.getAtomNameString() );

                break;

            case "C1'":

                an.setElement(ElementProvider.CARBON);
                an.setAtomNameString("C8'");
                logger.trace("Replacing Atom {} through {}.", "C1'",an.getAtomNameString() );

                break;

            case "O3'":

                an.setElement(ElementProvider.CARBON);
                an.setAtomNameString("C5'");
                logger.trace("Replacing Atom {} through {}.", "O3'",an.getAtomNameString() );

                break;

            case "P":

                an.setElement(ElementProvider.CARBON);
                an.setAtomNameString("C'");
                logger.trace("Replacing Atom {} through {}.", "P",an.getAtomNameString() );

                break;

                default:
                    break;
        }
    }



}
