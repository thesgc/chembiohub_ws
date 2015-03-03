CREATE TABLE compound_mols (
    molregno integer NOT NULL,
    ctab mol
);


ALTER TABLE ONLY compound_mols
    ADD CONSTRAINT compound_mols_pkey PRIMARY KEY (molregno);


CREATE INDEX rdkit_mol_idx ON compound_mols USING gist (ctab);

