# PNAtor

## Decription


## Quick start

- clone repro
- navigate to PNAGenerator.java in IDE of your choice
- choose your PDB Structure (DNA/RNA) 

```java
Structure structure = StructureParser.online()
        .pdbIdentifier("1OLD")
        .everything()
        .setOptions(StructureParserOptions.withSettings(
            StructureParserOptions.Setting.OMIT_HYDROGENS))
        .parse();

```
- execute the ```convertToPNAStructure(structure)``` method
- write reulst file
```java
try {
    StructureWriter.writeLeafSubstructureContainer(structure.getFirstModel(), 
        Paths.get("path/to/outputPDB.pdb"));
} catch (IOException e) {
    e.printStackTrace();
}
```
## Requirements
Make sure you have the Java 8 or later installed.

