create sequence molregno_id_seq;
alter table molecule_dictionary alter column molregno set default nextval('molregno_id_seq');
alter sequence molregno_id_seq owned by molecule_dictionary.molregno;