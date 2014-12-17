--
-- PostgreSQL database dump
--

SET statement_timeout = 0;
SET lock_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SET check_function_bodies = false;
SET client_min_messages = warning;

--
-- Name: plpgsql; Type: EXTENSION; Schema: -; Owner: 
--

CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;


--
-- Name: EXTENSION plpgsql; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION plpgsql IS 'PL/pgSQL procedural language';


--
-- Name: hstore; Type: EXTENSION; Schema: -; Owner: 
--

CREATE EXTENSION IF NOT EXISTS hstore WITH SCHEMA public;


--
-- Name: EXTENSION hstore; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION hstore IS 'data type for storing sets of (key, value) pairs';


--
-- Name: rdkit; Type: EXTENSION; Schema: -; Owner: 
--

CREATE EXTENSION IF NOT EXISTS rdkit WITH SCHEMA public;


--
-- Name: EXTENSION rdkit; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION rdkit IS 'Cheminformatics functionality for PostgreSQL.';


SET search_path = public, pg_catalog;

SET default_tablespace = '';

SET default_with_oids = false;

--
-- Name: action_type; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE action_type (
    action_type character varying(50) NOT NULL,
    description character varying(200) NOT NULL,
    parent_type character varying(50)
);


ALTER TABLE public.action_type OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN action_type.action_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN action_type.action_type IS 'Primary key. Type of action of the drug e.g., agonist, antagonist';


--
-- Name: COLUMN action_type.description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN action_type.description IS 'Description of how the action type is used';


--
-- Name: COLUMN action_type.parent_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN action_type.parent_type IS 'Higher-level grouping of action types e.g., positive vs negative action';


--
-- Name: activities; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE activities (
    activity_id bigint NOT NULL,
    assay_id integer NOT NULL,
    doc_id integer,
    record_id integer NOT NULL,
    molregno integer,
    standard_relation character varying(50),
    published_value numeric,
    published_units character varying(100),
    standard_value numeric,
    standard_units character varying(100),
    standard_flag smallint,
    standard_type character varying(250),
    activity_comment character varying(4000),
    published_type character varying(250),
    data_validity_comment character varying(30),
    potential_duplicate smallint,
    published_relation character varying(50),
    pchembl_value numeric(4,2),
    bao_endpoint character varying(11),
    uo_units character varying(10),
    qudt_units character varying(70),
    CONSTRAINT activities_potential_duplicate_check CHECK (((potential_duplicate = ANY (ARRAY[0, 1])) OR (potential_duplicate IS NULL))),
    CONSTRAINT activities_standard_flag_check CHECK (((standard_flag = ANY (ARRAY[0, 1])) OR (standard_flag IS NULL)))
);


ALTER TABLE public.activities OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN activities.activity_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.activity_id IS 'Unique ID for the activity row';


--
-- Name: COLUMN activities.assay_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.assay_id IS 'Foreign key to the assays table (containing the assay description)';


--
-- Name: COLUMN activities.doc_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.doc_id IS 'Foreign key to documents table (for quick lookup of publication details - can also link to documents through compound_records or assays table)';


--
-- Name: COLUMN activities.record_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.record_id IS 'Foreign key to the compound_records table (containing information on the compound tested)';


--
-- Name: COLUMN activities.molregno; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.molregno IS 'Foreign key to compounds table (for quick lookup of compound structure - can also link to compounds through compound_records table)';


--
-- Name: COLUMN activities.standard_relation; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.standard_relation IS 'Symbol constraining the activity value (e.g. >, <, =)';


--
-- Name: COLUMN activities.published_value; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.published_value IS 'Datapoint value as it appears in the original publication.';


--
-- Name: COLUMN activities.published_units; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.published_units IS 'Units of measurement as they appear in the original publication';


--
-- Name: COLUMN activities.standard_value; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.standard_value IS 'Same as PUBLISHED_VALUE but transformed to common units: e.g. mM concentrations converted to nM.';


--
-- Name: COLUMN activities.standard_units; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.standard_units IS 'Selected ''Standard'' units for data type: e.g. concentrations are in nM.';


--
-- Name: COLUMN activities.standard_flag; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.standard_flag IS 'Shows whether the standardised columns have been curated/set (1) or just default to the published data (0).';


--
-- Name: COLUMN activities.standard_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.standard_type IS 'Standardised version of the published_activity_type (e.g. IC50 rather than Ic-50/Ic50/ic50/ic-50)';


--
-- Name: COLUMN activities.activity_comment; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.activity_comment IS 'Describes non-numeric activities i.e. ''Slighty active'', ''Not determined''';


--
-- Name: COLUMN activities.published_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.published_type IS 'Type of end-point measurement: e.g. IC50, LD50, %%inhibition etc, as it appears in the original publication';


--
-- Name: COLUMN activities.data_validity_comment; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.data_validity_comment IS 'Comment reflecting whether the values for this activity measurement are likely to be correct - one of ''Manually validated'' (checked original paper and value is correct), ''Potential author error'' (value looks incorrect but is as reported in the original paper), ''Outside typical range'' (value seems too high/low to be correct e.g., negative IC50 value), ''Non standard unit type'' (units look incorrect for this activity type).';


--
-- Name: COLUMN activities.potential_duplicate; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.potential_duplicate IS 'Indicates whether the value is likely to be a repeat citation of a value reported in a previous ChEMBL paper, rather than a new, independent measurement.';


--
-- Name: COLUMN activities.published_relation; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.published_relation IS 'Symbol constraining the activity value (e.g. >, <, =), as it appears in the original publication';


--
-- Name: COLUMN activities.pchembl_value; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.pchembl_value IS 'Negative log of selected concentration-response activity values (IC50/EC50/XC50/AC50/Ki/Kd/Potency)';


--
-- Name: COLUMN activities.bao_endpoint; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.bao_endpoint IS 'ID for the corresponding result type in BioAssay Ontology (based on standard_type)';


--
-- Name: COLUMN activities.uo_units; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.uo_units IS 'ID for the corresponding unit in Unit Ontology (based on standard_units)';


--
-- Name: COLUMN activities.qudt_units; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activities.qudt_units IS 'ID for the corresponding unit in QUDT Ontology (based on standard_units)';


--
-- Name: activity_stds_lookup; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE activity_stds_lookup (
    std_act_id integer NOT NULL,
    standard_type character varying(250) NOT NULL,
    definition character varying(500),
    standard_units character varying(100) NOT NULL,
    normal_range_min numeric(24,12),
    normal_range_max numeric(24,12)
);


ALTER TABLE public.activity_stds_lookup OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN activity_stds_lookup.std_act_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activity_stds_lookup.std_act_id IS 'Primary key.';


--
-- Name: COLUMN activity_stds_lookup.standard_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activity_stds_lookup.standard_type IS 'The standard_type that other published_types in the activities table have been converted to.';


--
-- Name: COLUMN activity_stds_lookup.definition; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activity_stds_lookup.definition IS 'A description/definition of the standard_type.';


--
-- Name: COLUMN activity_stds_lookup.standard_units; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activity_stds_lookup.standard_units IS 'The units that are applied to this standard_type and to which other published_units are converted. Note a standard_type may have more than one allowable standard_unit and therefore multiple rows in this table.';


--
-- Name: COLUMN activity_stds_lookup.normal_range_min; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activity_stds_lookup.normal_range_min IS 'The lowest value for this activity type that is likely to be genuine. This is only an approximation, so lower genuine values may exist, but it may be desirable to validate these before using them. For a given standard_type/units, values in the activities table below this threshold are flagged with a data_validity_comment of ''Outside typical range''.';


--
-- Name: COLUMN activity_stds_lookup.normal_range_max; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN activity_stds_lookup.normal_range_max IS 'The highest value for this activity type that is likely to be genuine. This is only an approximation, so higher genuine values may exist, but it may be desirable to validate these before using them. For a given standard_type/units, values in the activities table above this threshold are flagged with a data_validity_comment of ''Outside typical range''.';


--
-- Name: assay_parameters; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE assay_parameters (
    assay_param_id integer NOT NULL,
    assay_id integer NOT NULL,
    parameter_type character varying(20) NOT NULL,
    parameter_value character varying(2000) NOT NULL,
    CONSTRAINT assay_parameters_assay_param_id_check CHECK ((assay_param_id >= 0))
);


ALTER TABLE public.assay_parameters OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN assay_parameters.assay_param_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assay_parameters.assay_param_id IS 'Numeric primary key';


--
-- Name: COLUMN assay_parameters.assay_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assay_parameters.assay_id IS 'Foreign key to assays table. The assay to which this parameter belongs';


--
-- Name: COLUMN assay_parameters.parameter_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assay_parameters.parameter_type IS 'Foreign key to parameter_type table, defining the meaning of the parameter';


--
-- Name: COLUMN assay_parameters.parameter_value; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assay_parameters.parameter_value IS 'The value of the particular parameter';


--
-- Name: assay_type; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE assay_type (
    assay_type character varying(1) NOT NULL,
    assay_desc character varying(250)
);


ALTER TABLE public.assay_type OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN assay_type.assay_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assay_type.assay_type IS 'Single character representing assay type';


--
-- Name: COLUMN assay_type.assay_desc; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assay_type.assay_desc IS 'Description of assay type';


--
-- Name: assays; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE assays (
    assay_id integer NOT NULL,
    doc_id integer NOT NULL,
    description character varying(4000),
    assay_type character varying(1),
    assay_test_type character varying(20),
    assay_category character varying(20),
    assay_organism character varying(250),
    assay_tax_id bigint,
    assay_strain character varying(200),
    assay_tissue character varying(100),
    assay_cell_type character varying(100),
    assay_subcellular_fraction character varying(100),
    tid integer,
    relationship_type character varying(1),
    confidence_score smallint,
    curated_by character varying(32),
    src_id smallint NOT NULL,
    src_assay_id character varying(50),
    chembl_id character varying(20) NOT NULL,
    cell_id integer,
    bao_format character varying(11),
    CONSTRAINT assays_assay_category_check CHECK (((assay_category)::text = ANY (ARRAY[('screening'::character varying)::text, ('panel'::character varying)::text, ('confirmatory'::character varying)::text, ('summary'::character varying)::text, ('other'::character varying)::text]))),
    CONSTRAINT assays_assay_tax_id_check CHECK ((assay_tax_id >= 0)),
    CONSTRAINT assays_assay_test_type_check CHECK (((assay_test_type)::text = ANY (ARRAY[('In vivo'::character varying)::text, ('In vitro'::character varying)::text, ('Ex vivo'::character varying)::text]))),
    CONSTRAINT assays_confidence_score_check CHECK ((confidence_score >= 0))
);


ALTER TABLE public.assays OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN assays.assay_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.assay_id IS 'Unique ID for the assay';


--
-- Name: COLUMN assays.doc_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.doc_id IS 'Foreign key to documents table';


--
-- Name: COLUMN assays.description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.description IS 'Description of the reported assay';


--
-- Name: COLUMN assays.assay_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.assay_type IS 'Assay classification, e.g. B=Binding assay, A=ADME assay, F=Functional assay';


--
-- Name: COLUMN assays.assay_test_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.assay_test_type IS 'Type of assay system (i.e., in vivo or in vitro)';


--
-- Name: COLUMN assays.assay_category; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.assay_category IS 'screening, confirmatory (ie: dose-response), summary, panel or other.';


--
-- Name: COLUMN assays.assay_organism; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.assay_organism IS 'Name of the organism for the assay system (e.g., the organism, tissue or cell line in which an assay was performed). May differ from the target organism (e.g., for a human protein expressed in non-human cells, or pathogen-infected human cells).';


--
-- Name: COLUMN assays.assay_tax_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.assay_tax_id IS 'NCBI tax ID for the assay organism.';


--
-- Name: COLUMN assays.assay_strain; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.assay_strain IS 'Name of specific strain of the assay organism used (where known)';


--
-- Name: COLUMN assays.assay_tissue; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.assay_tissue IS 'Name of tissue used in the assay system (e.g., for tissue-based assays) or from which the assay system was derived (e.g., for cell/subcellular fraction-based assays).';


--
-- Name: COLUMN assays.assay_cell_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.assay_cell_type IS 'Name of cell type or cell line used in the assay system (e.g., for cell-based assays).';


--
-- Name: COLUMN assays.assay_subcellular_fraction; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.assay_subcellular_fraction IS 'Name of subcellular fraction used in the assay system (e.g., microsomes, mitochondria).';


--
-- Name: COLUMN assays.tid; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.tid IS 'Target identifier to which this assay has been mapped. Foreign key to target_dictionary. FROM cbh_chembl_usr_15 onwards, an assay will have only a single target assigned.';


--
-- Name: COLUMN assays.relationship_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.relationship_type IS 'Flag indicating of the relationship between the reported target in the source document and the assigned target from TARGET_DICTIONARY. Foreign key to RELATIONSHIP_TYPE table.';


--
-- Name: COLUMN assays.confidence_score; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.confidence_score IS 'Confidence score, indicating how accurately the assigned target(s) represents the actually assay target. Foreign key to CONFIDENCE_SCORE table. 0 means uncurated/unassigned, 1 = low confidence to 9 = high confidence.';


--
-- Name: COLUMN assays.curated_by; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.curated_by IS 'Indicates the level of curation of the target assignment. Foreign key to curation_lookup table.';


--
-- Name: COLUMN assays.src_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.src_id IS 'Foreign key to source table';


--
-- Name: COLUMN assays.src_assay_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.src_assay_id IS 'Identifier for the assay in the source database/deposition (e.g., pubchem AID)';


--
-- Name: COLUMN assays.chembl_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.chembl_id IS 'ChEMBL identifier for this assay (for use on web interface etc)';


--
-- Name: COLUMN assays.cell_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.cell_id IS 'Foreign key to cell dictionary. The cell type or cell line used in the assay';


--
-- Name: COLUMN assays.bao_format; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN assays.bao_format IS 'ID for the corresponding format type in BioAssay Ontology (e.g., cell-based, biochemical, organism-based etc)';


--
-- Name: atc_classification; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE atc_classification (
    who_name character varying(150),
    level1 character varying(10),
    level2 character varying(10),
    level3 character varying(10),
    level4 character varying(10),
    level5 character varying(10) NOT NULL,
    who_id character varying(15),
    level1_description character varying(150),
    level2_description character varying(150),
    level3_description character varying(150),
    level4_description character varying(150)
);


ALTER TABLE public.atc_classification OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN atc_classification.who_name; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN atc_classification.who_name IS 'WHO/INN name for the compound';


--
-- Name: COLUMN atc_classification.level1; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN atc_classification.level1 IS 'First level of classification';


--
-- Name: COLUMN atc_classification.level2; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN atc_classification.level2 IS 'Second level of classification';


--
-- Name: COLUMN atc_classification.level3; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN atc_classification.level3 IS 'Third level of classification';


--
-- Name: COLUMN atc_classification.level4; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN atc_classification.level4 IS 'Fourth level of classification';


--
-- Name: COLUMN atc_classification.level5; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN atc_classification.level5 IS 'Complete ATC code for compound';


--
-- Name: COLUMN atc_classification.who_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN atc_classification.who_id IS 'WHO Identifier for compound';


--
-- Name: COLUMN atc_classification.level1_description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN atc_classification.level1_description IS 'Description of first level of classification';


--
-- Name: COLUMN atc_classification.level2_description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN atc_classification.level2_description IS 'Description of second level of classification';


--
-- Name: COLUMN atc_classification.level3_description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN atc_classification.level3_description IS 'Description of third level of classification';


--
-- Name: COLUMN atc_classification.level4_description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN atc_classification.level4_description IS 'Description of fourth level of classification';


--
-- Name: auth_group; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE auth_group (
    id integer NOT NULL,
    name character varying(80) NOT NULL
);


ALTER TABLE public.auth_group OWNER TO cbh_chembl_usr;

--
-- Name: auth_group_id_seq; Type: SEQUENCE; Schema: public; Owner: chembl
--

CREATE SEQUENCE auth_group_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.auth_group_id_seq OWNER TO cbh_chembl_usr;

--
-- Name: auth_group_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: chembl
--

ALTER SEQUENCE auth_group_id_seq OWNED BY auth_group.id;


--
-- Name: auth_group_permissions; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE auth_group_permissions (
    id integer NOT NULL,
    group_id integer NOT NULL,
    permission_id integer NOT NULL
);


ALTER TABLE public.auth_group_permissions OWNER TO cbh_chembl_usr;

--
-- Name: auth_group_permissions_id_seq; Type: SEQUENCE; Schema: public; Owner: chembl
--

CREATE SEQUENCE auth_group_permissions_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.auth_group_permissions_id_seq OWNER TO cbh_chembl_usr;

--
-- Name: auth_group_permissions_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: chembl
--

ALTER SEQUENCE auth_group_permissions_id_seq OWNED BY auth_group_permissions.id;


--
-- Name: auth_permission; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE auth_permission (
    id integer NOT NULL,
    name character varying(50) NOT NULL,
    content_type_id integer NOT NULL,
    codename character varying(100) NOT NULL
);


ALTER TABLE public.auth_permission OWNER TO cbh_chembl_usr;

--
-- Name: auth_permission_id_seq; Type: SEQUENCE; Schema: public; Owner: chembl
--

CREATE SEQUENCE auth_permission_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.auth_permission_id_seq OWNER TO cbh_chembl_usr;

--
-- Name: auth_permission_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: chembl
--

ALTER SEQUENCE auth_permission_id_seq OWNED BY auth_permission.id;


--
-- Name: auth_user; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE auth_user (
    id integer NOT NULL,
    password character varying(128) NOT NULL,
    last_login timestamp with time zone NOT NULL,
    is_superuser boolean NOT NULL,
    username character varying(30) NOT NULL,
    first_name character varying(30) NOT NULL,
    last_name character varying(30) NOT NULL,
    email character varying(75) NOT NULL,
    is_staff boolean NOT NULL,
    is_active boolean NOT NULL,
    date_joined timestamp with time zone NOT NULL
);


ALTER TABLE public.auth_user OWNER TO cbh_chembl_usr;

--
-- Name: auth_user_groups; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE auth_user_groups (
    id integer NOT NULL,
    user_id integer NOT NULL,
    group_id integer NOT NULL
);


ALTER TABLE public.auth_user_groups OWNER TO cbh_chembl_usr;

--
-- Name: auth_user_groups_id_seq; Type: SEQUENCE; Schema: public; Owner: chembl
--

CREATE SEQUENCE auth_user_groups_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.auth_user_groups_id_seq OWNER TO cbh_chembl_usr;

--
-- Name: auth_user_groups_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: chembl
--

ALTER SEQUENCE auth_user_groups_id_seq OWNED BY auth_user_groups.id;


--
-- Name: auth_user_id_seq; Type: SEQUENCE; Schema: public; Owner: chembl
--

CREATE SEQUENCE auth_user_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.auth_user_id_seq OWNER TO cbh_chembl_usr;

--
-- Name: auth_user_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: chembl
--

ALTER SEQUENCE auth_user_id_seq OWNED BY auth_user.id;


--
-- Name: auth_user_user_permissions; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE auth_user_user_permissions (
    id integer NOT NULL,
    user_id integer NOT NULL,
    permission_id integer NOT NULL
);


ALTER TABLE public.auth_user_user_permissions OWNER TO cbh_chembl_usr;

--
-- Name: auth_user_user_permissions_id_seq; Type: SEQUENCE; Schema: public; Owner: chembl
--

CREATE SEQUENCE auth_user_user_permissions_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.auth_user_user_permissions_id_seq OWNER TO cbh_chembl_usr;

--
-- Name: auth_user_user_permissions_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: chembl
--

ALTER SEQUENCE auth_user_user_permissions_id_seq OWNED BY auth_user_user_permissions.id;


--
-- Name: binding_sites; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE binding_sites (
    site_id integer NOT NULL,
    site_name character varying(200),
    tid integer
);


ALTER TABLE public.binding_sites OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN binding_sites.site_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN binding_sites.site_id IS 'Primary key. Unique identifier for a binding site in a given target.';


--
-- Name: COLUMN binding_sites.site_name; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN binding_sites.site_name IS 'Name/label for the binding site.';


--
-- Name: COLUMN binding_sites.tid; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN binding_sites.tid IS 'Foreign key to target_dictionary. Target on which the binding site is found.';


--
-- Name: bio_component_sequences; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE bio_component_sequences (
    component_id integer NOT NULL,
    component_type character varying(50) NOT NULL,
    description character varying(200),
    sequence text,
    sequence_md5sum character varying(32),
    tax_id bigint,
    organism character varying(150),
    CONSTRAINT bio_component_sequences_tax_id_check CHECK ((tax_id >= 0))
);


ALTER TABLE public.bio_component_sequences OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN bio_component_sequences.component_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN bio_component_sequences.component_id IS 'Primary key. Unique identifier for each of the molecular components of biotherapeutics in ChEMBL (e.g., antibody chains, recombinant proteins, synthetic peptides).';


--
-- Name: COLUMN bio_component_sequences.component_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN bio_component_sequences.component_type IS 'Type of molecular component (e.g., ''PROTEIN'',''DNA'',''RNA'').';


--
-- Name: COLUMN bio_component_sequences.description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN bio_component_sequences.description IS 'Description/name of molecular component.';


--
-- Name: COLUMN bio_component_sequences.sequence; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN bio_component_sequences.sequence IS 'Sequence of the biotherapeutic component.';


--
-- Name: COLUMN bio_component_sequences.sequence_md5sum; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN bio_component_sequences.sequence_md5sum IS 'MD5 checksum of the sequence.';


--
-- Name: COLUMN bio_component_sequences.tax_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN bio_component_sequences.tax_id IS 'NCBI tax ID for the species from which the sequence is derived. May be null for humanized monoclonal antibodies, synthetic peptides etc.';


--
-- Name: COLUMN bio_component_sequences.organism; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN bio_component_sequences.organism IS 'Name of the species from which the sequence is derived.';


--
-- Name: biotherapeutic_components; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE biotherapeutic_components (
    biocomp_id integer NOT NULL,
    molregno integer NOT NULL,
    component_id integer NOT NULL
);


ALTER TABLE public.biotherapeutic_components OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN biotherapeutic_components.biocomp_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN biotherapeutic_components.biocomp_id IS 'Primary key.';


--
-- Name: COLUMN biotherapeutic_components.molregno; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN biotherapeutic_components.molregno IS 'Foreign key to the biotherapeutics table, indicating which biotherapeutic the component is part of.';


--
-- Name: COLUMN biotherapeutic_components.component_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN biotherapeutic_components.component_id IS 'Foreign key to the bio_component_sequences table, indicating which component is part of the biotherapeutic.';


--
-- Name: biotherapeutics; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE biotherapeutics (
    molregno integer NOT NULL,
    description character varying(2000)
);


ALTER TABLE public.biotherapeutics OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN biotherapeutics.molregno; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN biotherapeutics.molregno IS 'Foreign key to molecule_dictionary';


--
-- Name: COLUMN biotherapeutics.description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN biotherapeutics.description IS 'Description of the biotherapeutic.';


--
-- Name: cbh_chembl_model_extension_cbhcompoundbatch; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE cbh_chembl_model_extension_cbhcompoundbatch (
    id integer NOT NULL,
    ctab text,
    editable_by hstore NOT NULL,
    viewable_by hstore NOT NULL,
    related_molregno_id integer,
    warnings hstore NOT NULL,
    custom_fields hstore NOT NULL
);


ALTER TABLE public.cbh_chembl_model_extension_cbhcompoundbatch OWNER TO cbh_chembl_usr;

--
-- Name: cbh_chembl_model_extension_cbhcompoundbatch_id_seq; Type: SEQUENCE; Schema: public; Owner: chembl
--

CREATE SEQUENCE cbh_chembl_model_extension_cbhcompoundbatch_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.cbh_chembl_model_extension_cbhcompoundbatch_id_seq OWNER TO cbh_chembl_usr;

--
-- Name: cbh_chembl_model_extension_cbhcompoundbatch_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: chembl
--

ALTER SEQUENCE cbh_chembl_model_extension_cbhcompoundbatch_id_seq OWNED BY cbh_chembl_model_extension_cbhcompoundbatch.id;


--
-- Name: cell_dictionary; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE cell_dictionary (
    cell_id integer NOT NULL,
    cell_name character varying(50) NOT NULL,
    cell_description character varying(200),
    cell_source_tissue character varying(50),
    cell_source_organism character varying(150),
    cell_source_tax_id bigint,
    clo_id character varying(11),
    efo_id character varying(12),
    cellosaurus_id character varying(15),
    CONSTRAINT cell_dictionary_cell_source_tax_id_check CHECK ((cell_source_tax_id >= 0))
);


ALTER TABLE public.cell_dictionary OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN cell_dictionary.cell_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN cell_dictionary.cell_id IS 'Primary key. Unique identifier for each cell line in the target_dictionary.';


--
-- Name: COLUMN cell_dictionary.cell_name; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN cell_dictionary.cell_name IS 'Name of each cell line (as used in the target_dicitonary pref_name).';


--
-- Name: COLUMN cell_dictionary.cell_description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN cell_dictionary.cell_description IS 'Longer description (where available) of the cell line.';


--
-- Name: COLUMN cell_dictionary.cell_source_tissue; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN cell_dictionary.cell_source_tissue IS 'Tissue from which the cell line is derived, where known.';


--
-- Name: COLUMN cell_dictionary.cell_source_organism; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN cell_dictionary.cell_source_organism IS 'Name of organism from which the cell line is derived.';


--
-- Name: COLUMN cell_dictionary.cell_source_tax_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN cell_dictionary.cell_source_tax_id IS 'NCBI tax ID of the organism from which the cell line is derived.';


--
-- Name: COLUMN cell_dictionary.clo_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN cell_dictionary.clo_id IS 'ID for the corresponding cell line in Cell Line Ontology';


--
-- Name: COLUMN cell_dictionary.efo_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN cell_dictionary.efo_id IS 'ID for the corresponding cell line in Experimental Factory Ontology';


--
-- Name: COLUMN cell_dictionary.cellosaurus_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN cell_dictionary.cellosaurus_id IS 'ID for the corresponding cell line in Cellosaurus Ontology';


--
-- Name: chembl_id_lookup; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE chembl_id_lookup (
    chembl_id character varying(20) NOT NULL,
    entity_type character varying(50),
    entity_id integer,
    status character varying(10) DEFAULT 'ACTIVE'::character varying,
    CONSTRAINT chembl_id_lookup_entity_type_check CHECK (((entity_type)::text = ANY (ARRAY[('ASSAY'::character varying)::text, ('COMPOUND'::character varying)::text, ('DOCUMENT'::character varying)::text, ('TARGET'::character varying)::text]))),
    CONSTRAINT chembl_id_lookup_status_check CHECK (((status)::text = ANY (ARRAY[('ACTIVE'::character varying)::text, ('INACTIVE'::character varying)::text, ('OBS'::character varying)::text])))
);


ALTER TABLE public.chembl_id_lookup OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN chembl_id_lookup.chembl_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN chembl_id_lookup.chembl_id IS 'ChEMBL identifier';


--
-- Name: COLUMN chembl_id_lookup.entity_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN chembl_id_lookup.entity_type IS 'Type of entity (e.g., COMPOUND, ASSAY, TARGET)';


--
-- Name: COLUMN chembl_id_lookup.entity_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN chembl_id_lookup.entity_id IS 'Primary key for that entity in corresponding table (e.g., molregno for compounds, tid for targets)';


--
-- Name: COLUMN chembl_id_lookup.status; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN chembl_id_lookup.status IS 'Indicates whether the status of the entity within the database - ACTIVE, INACTIVE (downgraded), OBS (obsolete/removed).';


--
-- Name: component_class; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE component_class (
    component_id integer NOT NULL,
    protein_class_id integer NOT NULL,
    comp_class_id integer NOT NULL
);


ALTER TABLE public.component_class OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN component_class.component_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_class.component_id IS 'Foreign key to component_sequences table.';


--
-- Name: COLUMN component_class.protein_class_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_class.protein_class_id IS 'Foreign key to the protein_classification table.';


--
-- Name: COLUMN component_class.comp_class_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_class.comp_class_id IS 'Primary key.';


--
-- Name: component_domains; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE component_domains (
    compd_id integer NOT NULL,
    domain_id integer,
    component_id integer NOT NULL,
    start_position integer,
    end_position integer,
    CONSTRAINT component_domains_end_position_check CHECK ((end_position >= 0)),
    CONSTRAINT component_domains_start_position_check CHECK ((start_position >= 0))
);


ALTER TABLE public.component_domains OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN component_domains.compd_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_domains.compd_id IS 'Primary key.';


--
-- Name: COLUMN component_domains.domain_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_domains.domain_id IS 'Foreign key to the domains table, indicating the domain that is contained in the associated molecular component.';


--
-- Name: COLUMN component_domains.component_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_domains.component_id IS 'Foreign key to the component_sequences table, indicating the molecular_component that has the given domain.';


--
-- Name: COLUMN component_domains.start_position; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_domains.start_position IS 'Start position of the domain within the sequence given in the component_sequences table.';


--
-- Name: COLUMN component_domains.end_position; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_domains.end_position IS 'End position of the domain within the sequence given in the component_sequences table.';


--
-- Name: component_sequences; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE component_sequences (
    component_id integer NOT NULL,
    component_type character varying(50),
    accession character varying(25),
    sequence text,
    sequence_md5sum character varying(32),
    description character varying(200),
    tax_id bigint,
    organism character varying(150),
    db_source character varying(25),
    db_version character varying(10),
    CONSTRAINT component_sequences_component_type_check CHECK (((component_type)::text = ANY (ARRAY[('PROTEIN'::character varying)::text, ('DNA'::character varying)::text, ('RNA'::character varying)::text]))),
    CONSTRAINT component_sequences_db_source_check CHECK (((db_source)::text = ANY (ARRAY[('Manual'::character varying)::text, ('SWISS-PROT'::character varying)::text, ('TREMBL'::character varying)::text]))),
    CONSTRAINT component_sequences_tax_id_check CHECK ((tax_id >= 0))
);


ALTER TABLE public.component_sequences OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN component_sequences.component_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_sequences.component_id IS 'Primary key. Unique identifier for the component.';


--
-- Name: COLUMN component_sequences.component_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_sequences.component_type IS 'Type of molecular component represented (e.g., ''PROTEIN'',''DNA'',''RNA'').';


--
-- Name: COLUMN component_sequences.accession; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_sequences.accession IS 'Accession for the sequence in the source database from which it was taken (e.g., UniProt accession for proteins).';


--
-- Name: COLUMN component_sequences.sequence; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_sequences.sequence IS 'A representative sequence for the molecular component, as given in the source sequence database (not necessarily the exact sequence used in the assay).';


--
-- Name: COLUMN component_sequences.sequence_md5sum; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_sequences.sequence_md5sum IS 'MD5 checksum of the sequence.';


--
-- Name: COLUMN component_sequences.description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_sequences.description IS 'Description/name for the molecular component, usually taken from the source sequence database.';


--
-- Name: COLUMN component_sequences.tax_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_sequences.tax_id IS 'NCBI tax ID for the sequence in the source database (i.e., species that the protein/nucleic acid sequence comes from).';


--
-- Name: COLUMN component_sequences.organism; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_sequences.organism IS 'Name of the organism the sequence comes from.';


--
-- Name: COLUMN component_sequences.db_source; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_sequences.db_source IS 'The name of the source sequence database from which sequences/accessions are taken. For UniProt proteins, this field indicates whether the sequence is from SWISS-PROT or TREMBL.';


--
-- Name: COLUMN component_sequences.db_version; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_sequences.db_version IS 'The version of the source sequence database from which sequences/accession were last updated.';


--
-- Name: component_synonyms; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE component_synonyms (
    compsyn_id integer NOT NULL,
    component_id integer NOT NULL,
    component_synonym character varying(500),
    syn_type character varying(20),
    CONSTRAINT component_synonyms_syn_type_check CHECK (((syn_type)::text = ANY (ARRAY[('HGNC_SYMBOL'::character varying)::text, ('GENE_SYMBOL'::character varying)::text, ('UNIPROT'::character varying)::text, ('MANUAL'::character varying)::text, ('OTHER'::character varying)::text, ('EC_NUMBER'::character varying)::text])))
);


ALTER TABLE public.component_synonyms OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN component_synonyms.compsyn_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_synonyms.compsyn_id IS 'Primary key.';


--
-- Name: COLUMN component_synonyms.component_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_synonyms.component_id IS 'Foreign key to the component_sequences table. The component to which this synonym applies.';


--
-- Name: COLUMN component_synonyms.component_synonym; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_synonyms.component_synonym IS 'The synonym for the component.';


--
-- Name: COLUMN component_synonyms.syn_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN component_synonyms.syn_type IS 'The type or origin of the synonym (e.g., GENE_SYMBOL).';


--
-- Name: compound_mols; Type: TABLE; Schema: public; Owner: vagrant; Tablespace: 
--

CREATE TABLE compound_mols (
    molregno integer NOT NULL,
    ctab mol
);


ALTER TABLE public.compound_mols OWNER TO vagrant;

--
-- Name: compound_properties; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE compound_properties (
    molregno integer NOT NULL,
    mw_freebase numeric(9,2),
    alogp numeric(9,2),
    hba smallint,
    hbd smallint,
    psa numeric(9,2),
    rtb smallint,
    ro3_pass character varying(3),
    num_ro5_violations smallint,
    med_chem_friendly character varying(3),
    acd_most_apka numeric(9,2),
    acd_most_bpka numeric(9,2),
    acd_logp numeric(9,2),
    acd_logd numeric(9,2),
    molecular_species character varying(50),
    full_mwt numeric(9,2),
    aromatic_rings smallint,
    heavy_atoms smallint,
    num_alerts smallint,
    qed_weighted numeric(3,2),
    mw_monoisotopic numeric(11,4),
    full_molformula character varying(100),
    CONSTRAINT compound_properties_aromatic_rings_check CHECK ((aromatic_rings >= 0)),
    CONSTRAINT compound_properties_full_mwt_check CHECK ((full_mwt >= (0)::numeric)),
    CONSTRAINT compound_properties_hba_check CHECK ((hba >= 0)),
    CONSTRAINT compound_properties_hbd_check CHECK ((hbd >= 0)),
    CONSTRAINT compound_properties_heavy_atoms_check CHECK ((heavy_atoms >= 0)),
    CONSTRAINT compound_properties_med_chem_friendly_check CHECK (((med_chem_friendly)::text = ANY (ARRAY[('Y'::character varying)::text, ('N'::character varying)::text]))),
    CONSTRAINT compound_properties_molecular_species_check CHECK (((molecular_species)::text = ANY (ARRAY[('ACID'::character varying)::text, ('BASE'::character varying)::text, ('ZWITTERION'::character varying)::text, ('NEUTRAL'::character varying)::text]))),
    CONSTRAINT compound_properties_mw_freebase_check CHECK ((mw_freebase >= (0)::numeric)),
    CONSTRAINT compound_properties_mw_monoisotopic_check CHECK ((mw_monoisotopic >= (0)::numeric)),
    CONSTRAINT compound_properties_num_alerts_check CHECK ((num_alerts >= 0)),
    CONSTRAINT compound_properties_num_ro5_violations_check CHECK (((num_ro5_violations >= 0) AND (num_ro5_violations = ANY (ARRAY[0, 1, 2, 3, 4])))),
    CONSTRAINT compound_properties_psa_check CHECK ((psa >= (0)::numeric)),
    CONSTRAINT compound_properties_qed_weighted_check CHECK ((qed_weighted >= (0)::numeric)),
    CONSTRAINT compound_properties_ro3_pass_check CHECK (((ro3_pass)::text = ANY (ARRAY[('Y'::character varying)::text, ('N'::character varying)::text]))),
    CONSTRAINT compound_properties_rtb_check CHECK ((rtb >= 0))
);


ALTER TABLE public.compound_properties OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN compound_properties.molregno; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.molregno IS 'Foreign key to compounds table (compound structure)';


--
-- Name: COLUMN compound_properties.mw_freebase; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.mw_freebase IS 'Molecular weight of parent compound';


--
-- Name: COLUMN compound_properties.alogp; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.alogp IS 'Calculated ALogP';


--
-- Name: COLUMN compound_properties.hba; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.hba IS 'Number hydrogen bond acceptors';


--
-- Name: COLUMN compound_properties.hbd; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.hbd IS 'Number hydrogen bond donors';


--
-- Name: COLUMN compound_properties.psa; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.psa IS 'Polar surface area';


--
-- Name: COLUMN compound_properties.rtb; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.rtb IS 'Number rotatable bonds';


--
-- Name: COLUMN compound_properties.ro3_pass; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.ro3_pass IS 'Indicates whether the compound passes the rule-of-three (mw < 300, logP < 3 etc)';


--
-- Name: COLUMN compound_properties.num_ro5_violations; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.num_ro5_violations IS 'Number of violations of rule-of-five';


--
-- Name: COLUMN compound_properties.med_chem_friendly; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.med_chem_friendly IS 'Indicates whether the compound is considered Med Chem friendly (Y/N)';


--
-- Name: COLUMN compound_properties.acd_most_apka; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.acd_most_apka IS 'The most acidic pKa calculated using ACDlabs v12.01';


--
-- Name: COLUMN compound_properties.acd_most_bpka; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.acd_most_bpka IS 'The most basic pKa calculated using ACDlabs v12.01';


--
-- Name: COLUMN compound_properties.acd_logp; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.acd_logp IS 'The calculated octanol/water partition coefficient using ACDlabs v12.01';


--
-- Name: COLUMN compound_properties.acd_logd; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.acd_logd IS 'The calculated octanol/water distribution coefficient at pH7.4 using ACDlabs v12.01';


--
-- Name: COLUMN compound_properties.molecular_species; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.molecular_species IS 'Indicates whether the compound is an acid/base/neutral';


--
-- Name: COLUMN compound_properties.full_mwt; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.full_mwt IS 'Molecular weight of the full compound including any salts';


--
-- Name: COLUMN compound_properties.aromatic_rings; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.aromatic_rings IS 'Number of aromatic rings';


--
-- Name: COLUMN compound_properties.heavy_atoms; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.heavy_atoms IS 'Number of heavy (non-hydrogen) atoms';


--
-- Name: COLUMN compound_properties.num_alerts; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.num_alerts IS 'Number of structural alerts (as defined by Brenk et al., ChemMedChem 2008)';


--
-- Name: COLUMN compound_properties.qed_weighted; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.qed_weighted IS 'Weighted quantitative estimate of drug likeness (as defined by Bickerton et al., Nature Chem 2012)';


--
-- Name: COLUMN compound_properties.mw_monoisotopic; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.mw_monoisotopic IS 'Monoisotopic parent molecular weight';


--
-- Name: COLUMN compound_properties.full_molformula; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_properties.full_molformula IS 'Molecular formula for the full compound (including any salt)';


--
-- Name: compound_records; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE compound_records (
    record_id integer NOT NULL,
    molregno integer,
    doc_id integer NOT NULL,
    compound_key character varying(250),
    compound_name character varying(4000),
    src_id smallint NOT NULL,
    src_compound_id character varying(150)
);


ALTER TABLE public.compound_records OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN compound_records.record_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_records.record_id IS 'Unique ID for a compound/record';


--
-- Name: COLUMN compound_records.molregno; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_records.molregno IS 'Foreign key to compounds table (compound structure)';


--
-- Name: COLUMN compound_records.doc_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_records.doc_id IS 'Foreign key to documents table';


--
-- Name: COLUMN compound_records.compound_key; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_records.compound_key IS 'Key text identifying this compound in the scientific document';


--
-- Name: COLUMN compound_records.compound_name; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_records.compound_name IS 'Name of this compound recorded in the scientific document';


--
-- Name: COLUMN compound_records.src_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_records.src_id IS 'Foreign key to source table';


--
-- Name: COLUMN compound_records.src_compound_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_records.src_compound_id IS 'Identifier for the compound in the source database (e.g., pubchem SID)';


--
-- Name: compound_structures; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE compound_structures (
    molregno integer NOT NULL,
    molfile text,
    standard_inchi character varying(4000),
    standard_inchi_key character varying(27) NOT NULL,
    canonical_smiles character varying(4000)
);


ALTER TABLE public.compound_structures OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN compound_structures.molregno; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_structures.molregno IS 'Internal Primary Key for the compound structure and foreign key to molecule_dictionary table';


--
-- Name: COLUMN compound_structures.molfile; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_structures.molfile IS 'MDL Connection table representation of compound';


--
-- Name: COLUMN compound_structures.standard_inchi; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_structures.standard_inchi IS 'IUPAC standard InChI for the compound';


--
-- Name: COLUMN compound_structures.standard_inchi_key; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_structures.standard_inchi_key IS 'IUPAC standard InChI key for the compound';


--
-- Name: COLUMN compound_structures.canonical_smiles; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN compound_structures.canonical_smiles IS 'Canonical smiles, generated using pipeline pilot';


--
-- Name: confidence_score_lookup; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE confidence_score_lookup (
    confidence_score smallint NOT NULL,
    description character varying(100) NOT NULL,
    target_mapping character varying(30) NOT NULL,
    CONSTRAINT confidence_score_lookup_confidence_score_check CHECK ((confidence_score >= 0))
);


ALTER TABLE public.confidence_score_lookup OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN confidence_score_lookup.confidence_score; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN confidence_score_lookup.confidence_score IS '0-9 score showing level of confidence in assignment of the precise molecular target of the assay';


--
-- Name: COLUMN confidence_score_lookup.description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN confidence_score_lookup.description IS 'Description of the target types assigned with each score';


--
-- Name: COLUMN confidence_score_lookup.target_mapping; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN confidence_score_lookup.target_mapping IS 'Short description of the target types assigned with each score';


--
-- Name: corsheaders_corsmodel; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE corsheaders_corsmodel (
    id integer NOT NULL,
    cors character varying(255) NOT NULL
);


ALTER TABLE public.corsheaders_corsmodel OWNER TO cbh_chembl_usr;

--
-- Name: corsheaders_corsmodel_id_seq; Type: SEQUENCE; Schema: public; Owner: chembl
--

CREATE SEQUENCE corsheaders_corsmodel_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.corsheaders_corsmodel_id_seq OWNER TO cbh_chembl_usr;

--
-- Name: corsheaders_corsmodel_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: chembl
--

ALTER SEQUENCE corsheaders_corsmodel_id_seq OWNED BY corsheaders_corsmodel.id;


--
-- Name: curation_lookup; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE curation_lookup (
    curated_by character varying(32) NOT NULL,
    description character varying(100) NOT NULL
);


ALTER TABLE public.curation_lookup OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN curation_lookup.curated_by; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN curation_lookup.curated_by IS 'Short description of the level of curation';


--
-- Name: COLUMN curation_lookup.description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN curation_lookup.description IS 'Definition of terms in the curated_by field.';


--
-- Name: data_validity_lookup; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE data_validity_lookup (
    data_validity_comment character varying(30) NOT NULL,
    description character varying(200)
);


ALTER TABLE public.data_validity_lookup OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN data_validity_lookup.data_validity_comment; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN data_validity_lookup.data_validity_comment IS 'Primary key. Short description of various types of errors/warnings applied to values in the activities table.';


--
-- Name: COLUMN data_validity_lookup.description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN data_validity_lookup.description IS 'Definition of the terms in the data_validity_comment field.';


--
-- Name: defined_daily_dose; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE defined_daily_dose (
    atc_code character varying(10) NOT NULL,
    ddd_value numeric(9,2),
    ddd_units character varying(20),
    ddd_admr character varying(30),
    ddd_comment character varying(400),
    ddd_id integer NOT NULL,
    CONSTRAINT defined_daily_dose_ddd_units_check CHECK (((ddd_units)::text = ANY (ARRAY[('LSU'::character varying)::text, ('MU'::character varying)::text, ('TU'::character varying)::text, ('U'::character varying)::text, ('g'::character varying)::text, ('mcg'::character varying)::text, ('mg'::character varying)::text, ('ml'::character varying)::text, ('mmol'::character varying)::text, ('tablet'::character varying)::text])))
);


ALTER TABLE public.defined_daily_dose OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN defined_daily_dose.atc_code; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN defined_daily_dose.atc_code IS 'ATC code for the compound (foreign key to ATC_CLASSIFICATION table)';


--
-- Name: COLUMN defined_daily_dose.ddd_value; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN defined_daily_dose.ddd_value IS 'Value of defined daily dose';


--
-- Name: COLUMN defined_daily_dose.ddd_units; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN defined_daily_dose.ddd_units IS 'Units of defined daily dose';


--
-- Name: COLUMN defined_daily_dose.ddd_admr; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN defined_daily_dose.ddd_admr IS 'Administration route for dose';


--
-- Name: COLUMN defined_daily_dose.ddd_comment; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN defined_daily_dose.ddd_comment IS 'Comment';


--
-- Name: COLUMN defined_daily_dose.ddd_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN defined_daily_dose.ddd_id IS 'Internal primary key';


--
-- Name: django_content_type; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE django_content_type (
    id integer NOT NULL,
    name character varying(100) NOT NULL,
    app_label character varying(100) NOT NULL,
    model character varying(100) NOT NULL
);


ALTER TABLE public.django_content_type OWNER TO cbh_chembl_usr;

--
-- Name: django_content_type_id_seq; Type: SEQUENCE; Schema: public; Owner: chembl
--

CREATE SEQUENCE django_content_type_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.django_content_type_id_seq OWNER TO cbh_chembl_usr;

--
-- Name: django_content_type_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: chembl
--

ALTER SEQUENCE django_content_type_id_seq OWNED BY django_content_type.id;


--
-- Name: django_session; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE django_session (
    session_key character varying(40) NOT NULL,
    session_data text NOT NULL,
    expire_date timestamp with time zone NOT NULL
);


ALTER TABLE public.django_session OWNER TO cbh_chembl_usr;

--
-- Name: django_site; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE django_site (
    id integer NOT NULL,
    domain character varying(100) NOT NULL,
    name character varying(50) NOT NULL
);


ALTER TABLE public.django_site OWNER TO cbh_chembl_usr;

--
-- Name: django_site_id_seq; Type: SEQUENCE; Schema: public; Owner: chembl
--

CREATE SEQUENCE django_site_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.django_site_id_seq OWNER TO cbh_chembl_usr;

--
-- Name: django_site_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: chembl
--

ALTER SEQUENCE django_site_id_seq OWNED BY django_site.id;


--
-- Name: docs; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE docs (
    doc_id integer NOT NULL,
    journal character varying(50),
    year smallint,
    volume character varying(50),
    issue character varying(50),
    first_page character varying(50),
    last_page character varying(50),
    pubmed_id bigint,
    doi character varying(50),
    chembl_id character varying(20) NOT NULL,
    title character varying(500),
    doc_type character varying(50) NOT NULL,
    authors character varying(4000),
    abstract text,
    CONSTRAINT docs_doc_type_check CHECK (((doc_type)::text = ANY (ARRAY[('PUBLICATION'::character varying)::text, ('BOOK'::character varying)::text, ('DATASET'::character varying)::text]))),
    CONSTRAINT docs_pubmed_id_check CHECK ((pubmed_id >= 0)),
    CONSTRAINT docs_year_check CHECK ((year >= 0))
);


ALTER TABLE public.docs OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN docs.doc_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN docs.doc_id IS 'Unique ID for the document';


--
-- Name: COLUMN docs.journal; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN docs.journal IS 'Abbreviated journal name for an article';


--
-- Name: COLUMN docs.year; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN docs.year IS 'Year of journal article publication';


--
-- Name: COLUMN docs.volume; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN docs.volume IS 'Volume of journal article';


--
-- Name: COLUMN docs.issue; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN docs.issue IS 'Issue of journal article';


--
-- Name: COLUMN docs.first_page; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN docs.first_page IS 'First page number of journal article';


--
-- Name: COLUMN docs.last_page; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN docs.last_page IS 'Last page number of journal article';


--
-- Name: COLUMN docs.pubmed_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN docs.pubmed_id IS 'NIH pubmed record ID, where available';


--
-- Name: COLUMN docs.doi; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN docs.doi IS 'Digital object identifier for this reference';


--
-- Name: COLUMN docs.chembl_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN docs.chembl_id IS 'ChEMBL identifier for this document (for use on web interface etc)';


--
-- Name: COLUMN docs.title; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN docs.title IS 'Document title (e.g., Publication title or description of dataset)';


--
-- Name: COLUMN docs.doc_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN docs.doc_type IS 'Type of the document (e.g., Publication, Deposited dataset)';


--
-- Name: COLUMN docs.authors; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN docs.authors IS 'For a deposited dataset, the authors carrying out the screening and/or submitting the dataset.';


--
-- Name: COLUMN docs.abstract; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN docs.abstract IS 'For a deposited dataset, a brief description of the dataset.';


--
-- Name: domains; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE domains (
    domain_id integer NOT NULL,
    domain_type character varying(20) NOT NULL,
    source_domain_id character varying(20) NOT NULL,
    domain_name character varying(20),
    domain_description character varying(500),
    CONSTRAINT domains_domain_type_check CHECK (((domain_type)::text = ANY (ARRAY[('Pfam-A'::character varying)::text, ('Pfam-B'::character varying)::text])))
);


ALTER TABLE public.domains OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN domains.domain_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN domains.domain_id IS 'Primary key. Unique identifier for each domain.';


--
-- Name: COLUMN domains.domain_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN domains.domain_type IS 'Indicates the source of the domain (e.g., Pfam).';


--
-- Name: COLUMN domains.source_domain_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN domains.source_domain_id IS 'Identifier for the domain in the source database (e.g., Pfam ID such as PF00001).';


--
-- Name: COLUMN domains.domain_name; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN domains.domain_name IS 'Name given to the domain in the source database (e.g., 7tm_1).';


--
-- Name: COLUMN domains.domain_description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN domains.domain_description IS 'Longer name or description for the domain.';


--
-- Name: drug_mechanism; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE drug_mechanism (
    mec_id integer NOT NULL,
    record_id integer NOT NULL,
    molregno integer,
    mechanism_of_action character varying(250),
    tid integer,
    site_id integer,
    action_type character varying(50),
    direct_interaction smallint,
    molecular_mechanism smallint,
    disease_efficacy smallint,
    mechanism_comment character varying(500),
    selectivity_comment character varying(100),
    binding_site_comment character varying(100),
    CONSTRAINT drug_mechanism_direct_interaction_check CHECK (((direct_interaction = ANY (ARRAY[0, 1])) OR (direct_interaction IS NULL))),
    CONSTRAINT drug_mechanism_disease_efficacy_check CHECK (((disease_efficacy = ANY (ARRAY[0, 1])) OR (disease_efficacy IS NULL))),
    CONSTRAINT drug_mechanism_molecular_mechanism_check CHECK (((molecular_mechanism = ANY (ARRAY[0, 1])) OR (molecular_mechanism IS NULL))),
    CONSTRAINT drug_mechanism_selectivity_comment_check CHECK (((selectivity_comment)::text = ANY (ARRAY[('Broad spectrum'::character varying)::text, ('EDG5 less relevant'::character varying)::text, ('M3 selective'::character varying)::text, ('Non-selective but type 5 receptor is overexpressed in Cushing''s disease'::character varying)::text, ('Selective'::character varying)::text, ('Selective for the brain omega-1 receptor (i.e. BZ1-type, i.e. alpha1/beta1/gamma2-GABA receptor)'::character varying)::text, ('Selectivity for types 2, 3 and 5'::character varying)::text, ('selectivity for beta-3 containing complexes'::character varying)::text])))
);


ALTER TABLE public.drug_mechanism OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN drug_mechanism.mec_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN drug_mechanism.mec_id IS 'Primary key for each drug mechanism of action';


--
-- Name: COLUMN drug_mechanism.record_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN drug_mechanism.record_id IS 'Record_id for the drug (foreign key to compound_records table)';


--
-- Name: COLUMN drug_mechanism.molregno; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN drug_mechanism.molregno IS 'Molregno for the drug (foreign key to molecule_dictionary table)';


--
-- Name: COLUMN drug_mechanism.mechanism_of_action; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN drug_mechanism.mechanism_of_action IS 'Description of the mechanism of action e.g., ''Phosphodiesterase 5 inhibitor''';


--
-- Name: COLUMN drug_mechanism.tid; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN drug_mechanism.tid IS 'Target associated with this mechanism of action (foreign key to target_dictionary table)';


--
-- Name: COLUMN drug_mechanism.site_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN drug_mechanism.site_id IS 'Binding site for the drug within the target (where known) - foreign key to binding_sites table';


--
-- Name: COLUMN drug_mechanism.action_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN drug_mechanism.action_type IS 'Type of action of the drug on the target e.g., agonist/antagonist etc (foreign key to action_type table)';


--
-- Name: COLUMN drug_mechanism.direct_interaction; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN drug_mechanism.direct_interaction IS 'Flag to show whether the molecule is believed to interact directly with the target (1 = yes, 0 = no)';


--
-- Name: COLUMN drug_mechanism.molecular_mechanism; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN drug_mechanism.molecular_mechanism IS 'Flag to show whether the mechanism of action describes the molecular target of the drug, rather than a higher-level physiological mechanism e.g., vasodilator (1 = yes, 0 = no)';


--
-- Name: COLUMN drug_mechanism.disease_efficacy; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN drug_mechanism.disease_efficacy IS 'Flag to show whether the target assigned is believed to play a role in the efficacy of the drug in the indication(s) for which it is approved (1 = yes, 0 = no)';


--
-- Name: COLUMN drug_mechanism.mechanism_comment; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN drug_mechanism.mechanism_comment IS 'Additional comments regarding the mechanism of action';


--
-- Name: COLUMN drug_mechanism.selectivity_comment; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN drug_mechanism.selectivity_comment IS 'Additional comments regarding the selectivity of the drug';


--
-- Name: COLUMN drug_mechanism.binding_site_comment; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN drug_mechanism.binding_site_comment IS 'Additional comments regarding the binding site of the drug';


--
-- Name: formulations; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE formulations (
    product_id character varying(30) NOT NULL,
    ingredient character varying(200),
    strength character varying(200),
    record_id integer NOT NULL,
    molregno integer,
    formulation_id integer NOT NULL
);


ALTER TABLE public.formulations OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN formulations.product_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN formulations.product_id IS 'Unique identifier of the product. FK to PRODUCTS';


--
-- Name: COLUMN formulations.ingredient; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN formulations.ingredient IS 'Name of the approved ingredient within the product';


--
-- Name: COLUMN formulations.strength; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN formulations.strength IS 'Dose strength';


--
-- Name: COLUMN formulations.record_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN formulations.record_id IS 'Foreign key to the compound_records table.';


--
-- Name: COLUMN formulations.molregno; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN formulations.molregno IS 'Unique identifier of the ingredient FK to MOLECULE_DICTIONARY';


--
-- Name: COLUMN formulations.formulation_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN formulations.formulation_id IS 'Primary key.';


--
-- Name: fps2_rdkit; Type: TABLE; Schema: public; Owner: vagrant; Tablespace: 
--

CREATE TABLE fps2_rdkit (
    molregno integer NOT NULL,
    torsionbv bfp,
    mfp2 bfp,
    ffp2 bfp,
    rdkfp bfp,
    atombv bfp,
    layeredfp bfp,
    maccsfp bfp
);


ALTER TABLE public.fps2_rdkit OWNER TO vagrant;

--
-- Name: fps_rdkit; Type: TABLE; Schema: public; Owner: vagrant; Tablespace: 
--

CREATE TABLE fps_rdkit (
    molregno integer NOT NULL,
    torsionbv bfp,
    mfp2 bfp,
    ffp2 bfp
);


ALTER TABLE public.fps_rdkit OWNER TO vagrant;

--
-- Name: ligand_eff; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE ligand_eff (
    activity_id bigint NOT NULL,
    bei numeric(9,2),
    sei numeric(9,2),
    le numeric(9,2),
    lle numeric(9,2),
    CONSTRAINT ligand_eff_bei_check CHECK ((bei >= (0)::numeric)),
    CONSTRAINT ligand_eff_le_check CHECK ((le >= (0)::numeric)),
    CONSTRAINT ligand_eff_sei_check CHECK ((sei >= (0)::numeric))
);


ALTER TABLE public.ligand_eff OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN ligand_eff.activity_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN ligand_eff.activity_id IS 'Link key to activities table';


--
-- Name: COLUMN ligand_eff.bei; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN ligand_eff.bei IS 'Binding Efficiency Index = p(XC50) *1000/MW_freebase';


--
-- Name: COLUMN ligand_eff.sei; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN ligand_eff.sei IS 'Surface Efficiency Index = p(XC50)*100/PSA';


--
-- Name: COLUMN ligand_eff.le; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN ligand_eff.le IS 'Ligand Efficiency = deltaG/heavy_atoms  [from the Hopkins DDT paper 2004]';


--
-- Name: COLUMN ligand_eff.lle; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN ligand_eff.lle IS 'Lipophilic Ligand Efficiency = -logKi-ALogP. [from Leeson NRDD 2007]';


--
-- Name: mechanism_refs; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE mechanism_refs (
    mecref_id integer NOT NULL,
    mec_id integer NOT NULL,
    ref_type character varying(50) NOT NULL,
    ref_id character varying(100),
    ref_url character varying(200),
    CONSTRAINT mechanism_refs_mecref_id_check CHECK ((mecref_id >= 0)),
    CONSTRAINT mechanism_refs_ref_type_check CHECK (((ref_type)::text = ANY (ARRAY[('ISBN'::character varying)::text, ('IUPHAR'::character varying)::text, ('DOI'::character varying)::text, ('EMA'::character varying)::text, ('PubMed'::character varying)::text, ('USPO'::character varying)::text, ('DailyMed'::character varying)::text, ('FDA'::character varying)::text, ('Expert'::character varying)::text, ('Other'::character varying)::text, ('InterPro'::character varying)::text, ('Wikipedia'::character varying)::text, ('UniProt'::character varying)::text, ('KEGG'::character varying)::text])))
);


ALTER TABLE public.mechanism_refs OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN mechanism_refs.mecref_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN mechanism_refs.mecref_id IS 'Primary key';


--
-- Name: COLUMN mechanism_refs.mec_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN mechanism_refs.mec_id IS 'Foreign key to drug_mechanism table - indicating the mechanism to which the references refer';


--
-- Name: COLUMN mechanism_refs.ref_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN mechanism_refs.ref_type IS 'Type/source of reference (e.g., ''PubMed'',''DailyMed'')';


--
-- Name: COLUMN mechanism_refs.ref_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN mechanism_refs.ref_id IS 'Identifier for the reference in the source (e.g., PubMed ID or DailyMed setid)';


--
-- Name: COLUMN mechanism_refs.ref_url; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN mechanism_refs.ref_url IS 'Full URL linking to the reference';


--
-- Name: molecule_atc_classification; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE molecule_atc_classification (
    mol_atc_id integer NOT NULL,
    level5 character varying(10) NOT NULL,
    molregno integer NOT NULL
);


ALTER TABLE public.molecule_atc_classification OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN molecule_atc_classification.mol_atc_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_atc_classification.mol_atc_id IS 'Primary key';


--
-- Name: COLUMN molecule_atc_classification.level5; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_atc_classification.level5 IS 'ATC code (foreign key to atc_classification table)';


--
-- Name: COLUMN molecule_atc_classification.molregno; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_atc_classification.molregno IS 'Drug to which the ATC code applies (foreign key to molecule_dictionary table)';


--
-- Name: molecule_dictionary; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE molecule_dictionary (
    molregno integer NOT NULL,
    pref_name character varying(255),
    chembl_id character varying(20) NOT NULL,
    max_phase smallint DEFAULT 0 NOT NULL,
    therapeutic_flag smallint DEFAULT 0 NOT NULL,
    dosed_ingredient smallint DEFAULT 0 NOT NULL,
    structure_type character varying(10) NOT NULL,
    chebi_par_id integer,
    molecule_type character varying(30),
    first_approval smallint,
    oral smallint DEFAULT 0 NOT NULL,
    parenteral smallint DEFAULT 0 NOT NULL,
    topical smallint DEFAULT 0 NOT NULL,
    black_box_warning smallint DEFAULT 0 NOT NULL,
    natural_product smallint DEFAULT (-1) NOT NULL,
    first_in_class smallint DEFAULT (-1) NOT NULL,
    chirality smallint DEFAULT (-1) NOT NULL,
    prodrug smallint DEFAULT (-1) NOT NULL,
    inorganic_flag smallint DEFAULT 0 NOT NULL,
    usan_year smallint,
    availability_type smallint,
    usan_stem character varying(50),
    polymer_flag smallint,
    usan_substem character varying(50),
    usan_stem_definition character varying(1000),
    indication_class character varying(1000),
    structure_key character varying(30),
    nomerge_reason character varying(2000),
    chebi_id integer,
    insert_date date,
    molfile_update date,
    downgrade_reason character varying(2000),
    checked_by character varying(2000),
    replacement_mrn integer,
    exclude integer,
    downgraded integer,
    nomerge integer,
    CONSTRAINT molecule_dictionary_availability_type_check CHECK ((availability_type = ANY (ARRAY[(-1), 0, 1, 2]))),
    CONSTRAINT molecule_dictionary_black_box_warning_check CHECK ((black_box_warning = ANY (ARRAY[0, 1, (-1)]))),
    CONSTRAINT molecule_dictionary_chebi_par_id_check CHECK ((chebi_par_id >= 0)),
    CONSTRAINT molecule_dictionary_chirality_check CHECK ((chirality = ANY (ARRAY[(-1), 0, 1, 2]))),
    CONSTRAINT molecule_dictionary_dosed_ingredient_check CHECK ((dosed_ingredient = ANY (ARRAY[0, 1]))),
    CONSTRAINT molecule_dictionary_first_approval_check CHECK ((first_approval >= 0)),
    CONSTRAINT molecule_dictionary_first_in_class_check CHECK ((first_in_class = ANY (ARRAY[0, 1, (-1)]))),
    CONSTRAINT molecule_dictionary_inorganic_flag_check CHECK ((inorganic_flag = ANY (ARRAY[0, 1, (-1)]))),
    CONSTRAINT molecule_dictionary_max_phase_check CHECK (((max_phase >= 0) AND (max_phase = ANY (ARRAY[0, 1, 2, 3, 4])))),
    CONSTRAINT molecule_dictionary_molecule_type_check CHECK (((molecule_type)::text = ANY (ARRAY[('Antibody'::character varying)::text, ('Cell'::character varying)::text, ('Enzyme'::character varying)::text, ('Oligonucleotide'::character varying)::text, ('Oligosaccharide'::character varying)::text, ('Protein'::character varying)::text, ('Small molecule'::character varying)::text, ('Unclassified'::character varying)::text, ('Unknown'::character varying)::text]))),
    CONSTRAINT molecule_dictionary_natural_product_check CHECK ((natural_product = ANY (ARRAY[0, 1, (-1)]))),
    CONSTRAINT molecule_dictionary_oral_check CHECK ((oral = ANY (ARRAY[0, 1]))),
    CONSTRAINT molecule_dictionary_parenteral_check CHECK ((parenteral = ANY (ARRAY[0, 1]))),
    CONSTRAINT molecule_dictionary_polymer_flag_check CHECK (((polymer_flag = ANY (ARRAY[0, 1])) OR (polymer_flag IS NULL))),
    CONSTRAINT molecule_dictionary_prodrug_check CHECK ((prodrug = ANY (ARRAY[0, 1, (-1)]))),
    CONSTRAINT molecule_dictionary_structure_type_check CHECK (((structure_type)::text = ANY (ARRAY[('NONE'::character varying)::text, ('MOL'::character varying)::text, ('SEQ'::character varying)::text, ('BOTH'::character varying)::text]))),
    CONSTRAINT molecule_dictionary_therapeutic_flag_check CHECK ((therapeutic_flag = ANY (ARRAY[0, 1]))),
    CONSTRAINT molecule_dictionary_topical_check CHECK ((topical = ANY (ARRAY[0, 1]))),
    CONSTRAINT molecule_dictionary_usan_year_check CHECK ((usan_year >= 0))
);


ALTER TABLE public.molecule_dictionary OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN molecule_dictionary.molregno; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.molregno IS 'Internal Primary Key for the molecule';


--
-- Name: COLUMN molecule_dictionary.pref_name; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.pref_name IS 'Preferred name for the molecule';


--
-- Name: COLUMN molecule_dictionary.chembl_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.chembl_id IS 'ChEMBL identifier for this compound (for use on web interface etc)';


--
-- Name: COLUMN molecule_dictionary.max_phase; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.max_phase IS 'Maximum phase of development reached for the compound (4 = approved). Null where max phase has not yet been assigned.';


--
-- Name: COLUMN molecule_dictionary.therapeutic_flag; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.therapeutic_flag IS 'Indicates that a drug has a therapeutic application (as opposed to e.g., an imaging agent, additive etc).';


--
-- Name: COLUMN molecule_dictionary.dosed_ingredient; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.dosed_ingredient IS 'Indicates that the drug is dosed in this form (e.g., a particular salt)';


--
-- Name: COLUMN molecule_dictionary.structure_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.structure_type IS 'Indications whether the molecule has a small molecule structure or a protein sequence (MOL indicates an entry in the compound_structures table, SEQ indications an entry in the protein_therapeutics table, NONE indicates an entry in neither table, e.g., structure unknown)';


--
-- Name: COLUMN molecule_dictionary.chebi_par_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.chebi_par_id IS 'Preferred ChEBI ID for the compound (where different from assigned)';


--
-- Name: COLUMN molecule_dictionary.molecule_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.molecule_type IS 'Type of molecule (Small molecule, Protein, Antibody, Oligosaccharide, Oligonucleotide, Cell, Unknown)';


--
-- Name: COLUMN molecule_dictionary.first_approval; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.first_approval IS 'Earliest known approval year for the molecule';


--
-- Name: COLUMN molecule_dictionary.oral; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.oral IS 'Indicates whether the drug is known to be administered orally.';


--
-- Name: COLUMN molecule_dictionary.parenteral; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.parenteral IS 'Indicates whether the drug is known to be administered parenterally';


--
-- Name: COLUMN molecule_dictionary.topical; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.topical IS 'Indicates whether the drug is known to be administered topically.';


--
-- Name: COLUMN molecule_dictionary.black_box_warning; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.black_box_warning IS 'Indicates that the drug has a black box warning';


--
-- Name: COLUMN molecule_dictionary.natural_product; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.natural_product IS 'Indicates whether the compound is natural product-derived (currently curated only for drugs)';


--
-- Name: COLUMN molecule_dictionary.first_in_class; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.first_in_class IS 'Indicates whether this is known to be the first compound of its class (e.g., acting on a particular target).';


--
-- Name: COLUMN molecule_dictionary.chirality; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.chirality IS 'Shows whether a drug is dosed as a racemic mixture (0), single stereoisomer (1) or is an achiral molecule (2)';


--
-- Name: COLUMN molecule_dictionary.prodrug; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.prodrug IS 'Indicates that the molecule is a pro-drug (see molecule hierarchy for active component, where known)';


--
-- Name: COLUMN molecule_dictionary.inorganic_flag; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.inorganic_flag IS 'Indicates whether the molecule is inorganic (i.e., containing only metal atoms and <2 carbon atoms)';


--
-- Name: COLUMN molecule_dictionary.usan_year; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.usan_year IS 'The year in which the application for a USAN/INN name was made';


--
-- Name: COLUMN molecule_dictionary.availability_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.availability_type IS 'The availability type for the drug (0 = discontinued, 1 = prescription only, 2 = over the counter)';


--
-- Name: COLUMN molecule_dictionary.usan_stem; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.usan_stem IS 'Where the compound has been assigned a USAN name, this indicates the stem, as described in the USAN_STEM table.';


--
-- Name: COLUMN molecule_dictionary.polymer_flag; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.polymer_flag IS 'Indicates whether a molecule is a small molecule polymer (e.g., polistyrex)';


--
-- Name: COLUMN molecule_dictionary.usan_substem; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.usan_substem IS 'Where the compound has been assigned a USAN name, this indicates the substem';


--
-- Name: COLUMN molecule_dictionary.usan_stem_definition; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.usan_stem_definition IS 'Definition of the USAN stem';


--
-- Name: COLUMN molecule_dictionary.indication_class; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_dictionary.indication_class IS 'Indication class(es) assigned to a drug in the USP dictionary';


--
-- Name: molecule_hierarchy; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE molecule_hierarchy (
    molregno integer NOT NULL,
    parent_molregno integer,
    active_molregno integer
);


ALTER TABLE public.molecule_hierarchy OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN molecule_hierarchy.molregno; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_hierarchy.molregno IS 'Foreign key to compounds table. This field holds a list of all of the ChEMBL compounds with associated data (e.g., activity information, approved drugs). Parent compounds that are generated only by removing salts, and which do not themselves have any associated data will not appear here.';


--
-- Name: COLUMN molecule_hierarchy.parent_molregno; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_hierarchy.parent_molregno IS 'Represents parent compound of molregno in first field (i.e., generated by removing salts). Where molregno and parent_molregno are same, the initial ChEMBL compound did not contain a salt component, or else could not be further processed for various reasons (e.g., inorganic mixture). Compounds which are only generated by removing salts will appear in this field only. Those which, themselves, have any associated data (e.g., activity data) or are launched drugs will also appear in the molregno field.';


--
-- Name: COLUMN molecule_hierarchy.active_molregno; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_hierarchy.active_molregno IS 'Where a compound is a pro-drug, this represents the active metabolite of the ''dosed'' compound given by parent_molregno. Where parent_molregno and active_molregno are the same, the compound is not currently known to be a pro-drug. ';


--
-- Name: molecule_synonyms; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE molecule_synonyms (
    molregno integer NOT NULL,
    syn_type character varying(50) NOT NULL,
    molsyn_id integer NOT NULL,
    res_stem_id integer,
    synonyms character varying(200)
);


ALTER TABLE public.molecule_synonyms OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN molecule_synonyms.molregno; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_synonyms.molregno IS 'Foreign key to molecule_dictionary';


--
-- Name: COLUMN molecule_synonyms.syn_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_synonyms.syn_type IS 'Type of name/synonym (e.g., TRADE_NAME, RESEARCH_CODE, USAN)';


--
-- Name: COLUMN molecule_synonyms.molsyn_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_synonyms.molsyn_id IS 'Primary key.';


--
-- Name: COLUMN molecule_synonyms.res_stem_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_synonyms.res_stem_id IS 'Foreign key to the research_stem table. Where a synonym is a research code, this links to further information about the company associated with that code.';


--
-- Name: COLUMN molecule_synonyms.synonyms; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN molecule_synonyms.synonyms IS 'Synonym for the compound';


--
-- Name: molregno_id_seq; Type: SEQUENCE; Schema: public; Owner: chembl
--

CREATE SEQUENCE molregno_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.molregno_id_seq OWNER TO cbh_chembl_usr;

--
-- Name: molregno_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: chembl
--

ALTER SEQUENCE molregno_id_seq OWNED BY molecule_dictionary.molregno;


--
-- Name: octmp_summary_id_seq; Type: SEQUENCE; Schema: public; Owner: chembl
--

CREATE SEQUENCE octmp_summary_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.octmp_summary_id_seq OWNER TO cbh_chembl_usr;

--
-- Name: octmp_summary; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE octmp_summary (
    id numeric(11,0) DEFAULT nextval('octmp_summary_id_seq'::regclass) NOT NULL,
    table_name character varying(50),
    table_created timestamp without time zone,
    query_md5 character varying(32)
);


ALTER TABLE public.octmp_summary OWNER TO cbh_chembl_usr;

--
-- Name: organism_class; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE organism_class (
    oc_id integer NOT NULL,
    tax_id bigint,
    l1 character varying(200),
    l2 character varying(200),
    l3 character varying(200),
    CONSTRAINT organism_class_tax_id_check CHECK ((tax_id >= 0))
);


ALTER TABLE public.organism_class OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN organism_class.oc_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN organism_class.oc_id IS 'Internal primary key';


--
-- Name: COLUMN organism_class.tax_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN organism_class.tax_id IS 'NCBI taxonomy ID for the organism (corresponding to tax_ids in assay2target and target_dictionary tables)';


--
-- Name: COLUMN organism_class.l1; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN organism_class.l1 IS 'Highest level classification (e.g., Eukaryotes, Bacteria, Fungi etc)';


--
-- Name: COLUMN organism_class.l2; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN organism_class.l2 IS 'Second level classification';


--
-- Name: COLUMN organism_class.l3; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN organism_class.l3 IS 'Third level classification';


--
-- Name: parameter_type; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE parameter_type (
    parameter_type character varying(20) NOT NULL,
    description character varying(2000)
);


ALTER TABLE public.parameter_type OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN parameter_type.parameter_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN parameter_type.parameter_type IS 'Short name for the type of parameter associated with an assay';


--
-- Name: COLUMN parameter_type.description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN parameter_type.description IS 'Description of the parameter type';


--
-- Name: predicted_binding_domains; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE predicted_binding_domains (
    predbind_id integer NOT NULL,
    activity_id bigint,
    site_id integer,
    prediction_method character varying(50),
    confidence character varying(10),
    CONSTRAINT predicted_binding_domains_confidence_check CHECK (((confidence)::text = ANY (ARRAY[('high'::character varying)::text, ('medium'::character varying)::text, ('low'::character varying)::text]))),
    CONSTRAINT predicted_binding_domains_prediction_method_check CHECK (((prediction_method)::text = ANY (ARRAY[('Manual'::character varying)::text, ('Multi domain'::character varying)::text, ('Single domain'::character varying)::text])))
);


ALTER TABLE public.predicted_binding_domains OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN predicted_binding_domains.predbind_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN predicted_binding_domains.predbind_id IS 'Primary key.';


--
-- Name: COLUMN predicted_binding_domains.activity_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN predicted_binding_domains.activity_id IS 'Foreign key to the activities table, indicating the compound/assay(+target) combination for which this prediction is made.';


--
-- Name: COLUMN predicted_binding_domains.site_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN predicted_binding_domains.site_id IS 'Foreign key to the binding_sites table, indicating the binding site (domain) that the compound is predicted to bind to.';


--
-- Name: COLUMN predicted_binding_domains.prediction_method; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN predicted_binding_domains.prediction_method IS 'The method used to assign the binding domain (e.g., ''Single domain'' where the protein has only 1 domain, ''Multi domain'' where the protein has multiple domains, but only 1 is known to bind small molecules in other proteins).';


--
-- Name: COLUMN predicted_binding_domains.confidence; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN predicted_binding_domains.confidence IS 'The level of confidence assigned to the prediction (high where the protein has only 1 domain, medium where the compound has multiple domains, but only 1 known small molecule-binding domain).';


--
-- Name: products; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE products (
    dosage_form character varying(200),
    route character varying(200),
    trade_name character varying(200),
    approval_date date,
    ad_type character varying(5),
    oral smallint,
    topical smallint,
    parenteral smallint,
    black_box_warning smallint,
    applicant_full_name character varying(200),
    innovator_company smallint,
    product_id character varying(30) NOT NULL,
    nda_type character varying(10),
    CONSTRAINT products_ad_type_check CHECK (((ad_type)::text = ANY (ARRAY[('OTC'::character varying)::text, ('RX'::character varying)::text, ('DISCN'::character varying)::text]))),
    CONSTRAINT products_black_box_warning_check CHECK (((black_box_warning = ANY (ARRAY[0, 1])) OR (black_box_warning IS NULL))),
    CONSTRAINT products_innovator_company_check CHECK (((innovator_company = ANY (ARRAY[0, 1])) OR (innovator_company IS NULL))),
    CONSTRAINT products_nda_type_check CHECK (((nda_type)::text = ANY (ARRAY[('A'::character varying)::text, ('N'::character varying)::text]))),
    CONSTRAINT products_oral_check CHECK (((oral = ANY (ARRAY[0, 1])) OR (oral IS NULL))),
    CONSTRAINT products_parenteral_check CHECK (((parenteral = ANY (ARRAY[0, 1])) OR (parenteral IS NULL))),
    CONSTRAINT products_topical_check CHECK (((topical = ANY (ARRAY[0, 1])) OR (topical IS NULL)))
);


ALTER TABLE public.products OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN products.dosage_form; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN products.dosage_form IS 'The dosage form of the product (e.g., tablet, capsule etc)';


--
-- Name: COLUMN products.route; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN products.route IS 'The administration route of the product (e.g., oral, injection etc)';


--
-- Name: COLUMN products.trade_name; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN products.trade_name IS 'The trade name for the product';


--
-- Name: COLUMN products.approval_date; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN products.approval_date IS 'The FDA approval date for the product (not necessarily first approval of the active ingredient)';


--
-- Name: COLUMN products.ad_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN products.ad_type IS 'RX = prescription, OTC = over the counter, DISCN = discontinued';


--
-- Name: COLUMN products.oral; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN products.oral IS 'Flag to show whether product is orally delivered';


--
-- Name: COLUMN products.topical; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN products.topical IS 'Flag to show whether product is topically delivered';


--
-- Name: COLUMN products.parenteral; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN products.parenteral IS 'Flag to show whether product is parenterally delivered';


--
-- Name: COLUMN products.black_box_warning; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN products.black_box_warning IS 'Flag to show whether the product label has a black box warning';


--
-- Name: COLUMN products.applicant_full_name; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN products.applicant_full_name IS 'Name of the company applying for FDA approval';


--
-- Name: COLUMN products.innovator_company; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN products.innovator_company IS 'Flag to show whether the applicant is the innovator of the product';


--
-- Name: COLUMN products.product_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN products.product_id IS 'FDA application number for the product';


--
-- Name: COLUMN products.nda_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN products.nda_type IS 'New Drug Application Type. The type of new drug application approval.  New Drug Applications (NDA or innovator)  are "N".   Abbreviated New Drug Applications (ANDA or generic) are "A".';


--
-- Name: protein_class_synonyms; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE protein_class_synonyms (
    protclasssyn_id integer NOT NULL,
    protein_class_id integer NOT NULL,
    protein_class_synonym character varying(1000),
    syn_type character varying(20),
    CONSTRAINT protein_class_synonyms_syn_type_check CHECK (((syn_type)::text = ANY (ARRAY[('CHEMBL'::character varying)::text, ('CONCEPT_WIKI'::character varying)::text, ('UMLS'::character varying)::text, ('CW_XREF'::character varying)::text, ('MESH_XREF'::character varying)::text])))
);


ALTER TABLE public.protein_class_synonyms OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN protein_class_synonyms.protclasssyn_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_class_synonyms.protclasssyn_id IS 'Primary key.';


--
-- Name: COLUMN protein_class_synonyms.protein_class_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_class_synonyms.protein_class_id IS 'Foreign key to the PROTEIN_CLASSIFICATION table. The protein_class to which this synonym applies.';


--
-- Name: COLUMN protein_class_synonyms.protein_class_synonym; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_class_synonyms.protein_class_synonym IS 'The synonym for the protein class.';


--
-- Name: COLUMN protein_class_synonyms.syn_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_class_synonyms.syn_type IS 'The type or origin of the synonym (e.g., ChEMBL, Concept Wiki, UMLS).';


--
-- Name: protein_classification; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE protein_classification (
    protein_class_id integer NOT NULL,
    parent_id integer,
    pref_name character varying(500),
    short_name character varying(50),
    protein_class_desc character varying(410) NOT NULL,
    definition character varying(4000),
    class_level integer,
    CONSTRAINT protein_classification_class_level_check CHECK (((class_level >= 0) AND (class_level = ANY (ARRAY[0, 1, 2, 3, 4, 5, 6, 7, 8])))),
    CONSTRAINT protein_classification_parent_id_check CHECK ((parent_id >= 0))
);


ALTER TABLE public.protein_classification OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN protein_classification.protein_class_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_classification.protein_class_id IS 'Primary key. Unique identifier for each protein family classification.';


--
-- Name: COLUMN protein_classification.parent_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_classification.parent_id IS 'Protein_class_id for the parent of this protein family.';


--
-- Name: COLUMN protein_classification.pref_name; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_classification.pref_name IS 'Preferred/full name for this protein family.';


--
-- Name: COLUMN protein_classification.short_name; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_classification.short_name IS 'Short/abbreviated name for this protein family (not necessarily unique).';


--
-- Name: COLUMN protein_classification.protein_class_desc; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_classification.protein_class_desc IS 'Concatenated description of each classification for searching purposes etc.';


--
-- Name: COLUMN protein_classification.definition; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_classification.definition IS 'Definition of the protein family.';


--
-- Name: protein_family_classification; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE protein_family_classification (
    protein_class_id integer NOT NULL,
    protein_class_desc character varying(810) NOT NULL,
    l1 character varying(100) NOT NULL,
    l2 character varying(100),
    l3 character varying(100),
    l4 character varying(100),
    l5 character varying(100),
    l6 character varying(100),
    l7 character varying(100),
    l8 character varying(100)
);


ALTER TABLE public.protein_family_classification OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN protein_family_classification.protein_class_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_family_classification.protein_class_id IS 'Primary key. Unique identifier for each classification.';


--
-- Name: COLUMN protein_family_classification.protein_class_desc; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_family_classification.protein_class_desc IS 'Concatenated description of each classification for searching purposes etc.';


--
-- Name: COLUMN protein_family_classification.l1; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_family_classification.l1 IS 'First level classification (e.g., Enzyme, Transporter, Ion Channel).';


--
-- Name: COLUMN protein_family_classification.l2; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_family_classification.l2 IS 'Second level classification.';


--
-- Name: COLUMN protein_family_classification.l3; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_family_classification.l3 IS 'Third level classification.';


--
-- Name: COLUMN protein_family_classification.l4; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_family_classification.l4 IS 'Fourth level classification.';


--
-- Name: COLUMN protein_family_classification.l5; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_family_classification.l5 IS 'Fifth level classification.';


--
-- Name: COLUMN protein_family_classification.l6; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_family_classification.l6 IS 'Sixth level classification.';


--
-- Name: COLUMN protein_family_classification.l7; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_family_classification.l7 IS 'Seventh level classification.';


--
-- Name: COLUMN protein_family_classification.l8; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN protein_family_classification.l8 IS 'Eighth level classification.';


--
-- Name: relationship_type; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE relationship_type (
    relationship_type character varying(1) NOT NULL,
    relationship_desc character varying(250)
);


ALTER TABLE public.relationship_type OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN relationship_type.relationship_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN relationship_type.relationship_type IS 'Relationship_type flag used in the assay2target table';


--
-- Name: COLUMN relationship_type.relationship_desc; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN relationship_type.relationship_desc IS 'Description of relationship_type flags';


--
-- Name: research_companies; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE research_companies (
    co_stem_id integer NOT NULL,
    res_stem_id integer,
    company character varying(100),
    country character varying(50),
    previous_company character varying(100)
);


ALTER TABLE public.research_companies OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN research_companies.co_stem_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN research_companies.co_stem_id IS 'Primary key.';


--
-- Name: COLUMN research_companies.res_stem_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN research_companies.res_stem_id IS 'Foreign key to research_stem table.';


--
-- Name: COLUMN research_companies.company; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN research_companies.company IS 'Name of current company associated with this research code stem.';


--
-- Name: COLUMN research_companies.country; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN research_companies.country IS 'Country in which the company uses this research code stem.';


--
-- Name: COLUMN research_companies.previous_company; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN research_companies.previous_company IS 'Previous name of the company associated with this research code stem (e.g., if the company has undergone acquisitions/mergers).';


--
-- Name: research_stem; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE research_stem (
    res_stem_id integer NOT NULL,
    research_stem character varying(20)
);


ALTER TABLE public.research_stem OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN research_stem.res_stem_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN research_stem.res_stem_id IS 'Primary key. Unique ID for each research code stem.';


--
-- Name: COLUMN research_stem.research_stem; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN research_stem.research_stem IS 'The actual stem/prefix used in the research code.';


--
-- Name: site_components; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE site_components (
    sitecomp_id integer NOT NULL,
    site_id integer NOT NULL,
    component_id integer,
    domain_id integer,
    site_residues character varying(2000)
);


ALTER TABLE public.site_components OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN site_components.sitecomp_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN site_components.sitecomp_id IS 'Primary key.';


--
-- Name: COLUMN site_components.site_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN site_components.site_id IS 'Foreign key to binding_sites table.';


--
-- Name: COLUMN site_components.component_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN site_components.component_id IS 'Foreign key to the component_sequences table, indicating which molecular component of the target is involved in the binding site.';


--
-- Name: COLUMN site_components.domain_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN site_components.domain_id IS 'Foreign key to the domains table, indicating which domain of the given molecular component is involved in the binding site (where not known, the domain_id may be null).';


--
-- Name: COLUMN site_components.site_residues; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN site_components.site_residues IS 'List of residues from the given molecular component that make up the binding site (where not know, will be null).';


--
-- Name: source; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE source (
    src_id smallint NOT NULL,
    src_description character varying(500),
    src_short_name character varying(20)
);


ALTER TABLE public.source OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN source.src_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN source.src_id IS 'Identifier for each source (used in compound_records and assays tables)';


--
-- Name: COLUMN source.src_description; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN source.src_description IS 'Description of the data source';


--
-- Name: COLUMN source.src_short_name; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN source.src_short_name IS 'A short name for each data source, for display purposes';


--
-- Name: south_migrationhistory; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE south_migrationhistory (
    id integer NOT NULL,
    app_name character varying(255) NOT NULL,
    migration character varying(255) NOT NULL,
    applied timestamp with time zone NOT NULL
);


ALTER TABLE public.south_migrationhistory OWNER TO cbh_chembl_usr;

--
-- Name: south_migrationhistory_id_seq; Type: SEQUENCE; Schema: public; Owner: chembl
--

CREATE SEQUENCE south_migrationhistory_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.south_migrationhistory_id_seq OWNER TO cbh_chembl_usr;

--
-- Name: south_migrationhistory_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: chembl
--

ALTER SEQUENCE south_migrationhistory_id_seq OWNED BY south_migrationhistory.id;


--
-- Name: target_components; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE target_components (
    tid integer NOT NULL,
    component_id integer NOT NULL,
    targcomp_id integer NOT NULL,
    homologue smallint DEFAULT 0 NOT NULL,
    CONSTRAINT target_components_homologue_check CHECK (((homologue >= 0) AND (homologue = ANY (ARRAY[0, 1, 2]))))
);


ALTER TABLE public.target_components OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN target_components.tid; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_components.tid IS 'Foreign key to the target_dictionary, indicating the target to which the components belong.';


--
-- Name: COLUMN target_components.component_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_components.component_id IS 'Foreign key to the component_sequences table, indicating which components belong to the target.';


--
-- Name: COLUMN target_components.targcomp_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_components.targcomp_id IS 'Primary key.';


--
-- Name: COLUMN target_components.homologue; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_components.homologue IS 'Indicates that the given component is a homologue of the correct component (e.g., from a different species) when set to 1. This may be the case if the sequence for the correct protein/nucleic acid cannot be found in sequence databases.';


--
-- Name: target_dictionary; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE target_dictionary (
    tid integer NOT NULL,
    target_type character varying(30),
    pref_name character varying(200),
    tax_id bigint,
    organism character varying(150),
    chembl_id character varying(20) NOT NULL,
    species_group_flag smallint,
    CONSTRAINT target_dictionary_species_group_flag_check CHECK (((species_group_flag = ANY (ARRAY[0, 1])) OR (species_group_flag IS NULL))),
    CONSTRAINT target_dictionary_tax_id_check CHECK ((tax_id >= 0))
);


ALTER TABLE public.target_dictionary OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN target_dictionary.tid; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_dictionary.tid IS 'Unique ID for the target';


--
-- Name: COLUMN target_dictionary.target_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_dictionary.target_type IS 'Describes whether target is a protein, an organism, a tissue etc. Foreign key to TARGET_TYPE table.';


--
-- Name: COLUMN target_dictionary.pref_name; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_dictionary.pref_name IS 'Preferred target name: manually curated';


--
-- Name: COLUMN target_dictionary.tax_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_dictionary.tax_id IS 'NCBI taxonomy id of target';


--
-- Name: COLUMN target_dictionary.organism; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_dictionary.organism IS 'Source organism of molecuar target or tissue, or the target organism if compound activity is reported in an organism rather than a protein or tissue';


--
-- Name: COLUMN target_dictionary.chembl_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_dictionary.chembl_id IS 'ChEMBL identifier for this target (for use on web interface etc)';


--
-- Name: COLUMN target_dictionary.species_group_flag; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_dictionary.species_group_flag IS 'Flag to indicate whether the target represents a group of species, rather than an individual species (e.g., ''Bacterial DHFR''). Where set to 1, indicates that any associated target components will be a representative, rather than a comprehensive set.';


--
-- Name: target_relations; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE target_relations (
    tid integer NOT NULL,
    relationship character varying(20) NOT NULL,
    related_tid integer NOT NULL,
    targrel_id integer NOT NULL,
    CONSTRAINT target_relations_relationship_check CHECK (((relationship)::text = ANY (ARRAY[('EQUIVALENT TO'::character varying)::text, ('OVERLAPS WITH'::character varying)::text, ('SUBSET OF'::character varying)::text, ('SUPERSET OF'::character varying)::text]))),
    CONSTRAINT target_relations_targrel_id_check CHECK ((targrel_id >= 0))
);


ALTER TABLE public.target_relations OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN target_relations.tid; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_relations.tid IS 'Identifier for target of interest (foreign key to target_dictionary table)';


--
-- Name: COLUMN target_relations.relationship; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_relations.relationship IS 'Relationship between two targets (e.g., SUBSET OF, SUPERSET OF, OVERLAPS WITH)';


--
-- Name: COLUMN target_relations.related_tid; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_relations.related_tid IS 'Identifier for the target that is related to the target of interest (foreign key to target_dicitionary table)';


--
-- Name: COLUMN target_relations.targrel_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_relations.targrel_id IS 'Primary key';


--
-- Name: target_type; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE target_type (
    target_type character varying(30) NOT NULL,
    target_desc character varying(250),
    parent_type character varying(25)
);


ALTER TABLE public.target_type OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN target_type.target_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_type.target_type IS 'Target type (as used in target dictionary)';


--
-- Name: COLUMN target_type.target_desc; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_type.target_desc IS 'Description of target type';


--
-- Name: COLUMN target_type.parent_type; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN target_type.parent_type IS 'Higher level classification of target_type, allowing grouping of e.g., all ''PROTEIN'' targets, all ''NON-MOLECULAR'' targets etc.';


--
-- Name: tastypie_apiaccess; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE tastypie_apiaccess (
    id integer NOT NULL,
    identifier character varying(255) NOT NULL,
    url character varying(255) NOT NULL,
    request_method character varying(10) NOT NULL,
    accessed integer NOT NULL,
    CONSTRAINT tastypie_apiaccess_accessed_check CHECK ((accessed >= 0))
);


ALTER TABLE public.tastypie_apiaccess OWNER TO cbh_chembl_usr;

--
-- Name: tastypie_apiaccess_id_seq; Type: SEQUENCE; Schema: public; Owner: chembl
--

CREATE SEQUENCE tastypie_apiaccess_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.tastypie_apiaccess_id_seq OWNER TO cbh_chembl_usr;

--
-- Name: tastypie_apiaccess_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: chembl
--

ALTER SEQUENCE tastypie_apiaccess_id_seq OWNED BY tastypie_apiaccess.id;


--
-- Name: usan_stems; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE usan_stems (
    usan_stem_id integer NOT NULL,
    stem character varying(100) NOT NULL,
    subgroup character varying(100) NOT NULL,
    annotation character varying(2000),
    stem_class character varying(100),
    major_class character varying(100),
    who_extra smallint DEFAULT 0,
    CONSTRAINT usan_stems_major_class_check CHECK (((major_class)::text = ANY (ARRAY[('GPCR'::character varying)::text, ('NR'::character varying)::text, ('PDE'::character varying)::text, ('ion channel'::character varying)::text, ('kinase'::character varying)::text, ('protease'::character varying)::text]))),
    CONSTRAINT usan_stems_stem_class_check CHECK (((stem_class)::text = ANY (ARRAY[('Suffix'::character varying)::text, ('Prefix'::character varying)::text, ('Infix'::character varying)::text]))),
    CONSTRAINT usan_stems_usan_stem_id_check CHECK ((usan_stem_id >= 0)),
    CONSTRAINT usan_stems_who_extra_check CHECK (((who_extra = ANY (ARRAY[0, 1])) OR (who_extra IS NULL)))
);


ALTER TABLE public.usan_stems OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN usan_stems.usan_stem_id; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN usan_stems.usan_stem_id IS 'Numeric primary key.';


--
-- Name: COLUMN usan_stems.stem; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN usan_stems.stem IS 'Stem defined for use in United States Adopted Names.';


--
-- Name: COLUMN usan_stems.subgroup; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN usan_stems.subgroup IS 'More specific subgroup of the stem defined for use in United States Adopted Names.';


--
-- Name: COLUMN usan_stems.annotation; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN usan_stems.annotation IS 'Meaning of the stem (e.g., the class of compound it applies to).';


--
-- Name: COLUMN usan_stems.stem_class; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN usan_stems.stem_class IS 'Indicates whether stem is used as a Prefix/Infix/Suffix.';


--
-- Name: COLUMN usan_stems.major_class; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN usan_stems.major_class IS 'Protein family targeted by compounds of this class (e.g., GPCR/Ion channel/Protease) where known/applicable.';


--
-- Name: COLUMN usan_stems.who_extra; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN usan_stems.who_extra IS 'Stem not represented in USAN list, but added from WHO INN stem list (where set to 1).';


--
-- Name: version; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE version (
    name character varying(20) NOT NULL,
    creation_date date,
    comments character varying(2000)
);


ALTER TABLE public.version OWNER TO cbh_chembl_usr;

--
-- Name: COLUMN version.name; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN version.name IS 'Name of release version';


--
-- Name: COLUMN version.creation_date; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN version.creation_date IS 'Date database created';


--
-- Name: COLUMN version.comments; Type: COMMENT; Schema: public; Owner: chembl
--

COMMENT ON COLUMN version.comments IS 'Description of release version';


--
-- Name: ws_cache; Type: TABLE; Schema: public; Owner: chembl; Tablespace: 
--

CREATE TABLE ws_cache (
    cache_key character varying(255) NOT NULL,
    value text NOT NULL,
    expires timestamp with time zone NOT NULL
);


ALTER TABLE public.ws_cache OWNER TO cbh_chembl_usr;

--
-- Name: id; Type: DEFAULT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY auth_group ALTER COLUMN id SET DEFAULT nextval('auth_group_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY auth_group_permissions ALTER COLUMN id SET DEFAULT nextval('auth_group_permissions_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY auth_permission ALTER COLUMN id SET DEFAULT nextval('auth_permission_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY auth_user ALTER COLUMN id SET DEFAULT nextval('auth_user_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY auth_user_groups ALTER COLUMN id SET DEFAULT nextval('auth_user_groups_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY auth_user_user_permissions ALTER COLUMN id SET DEFAULT nextval('auth_user_user_permissions_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY cbh_chembl_model_extension_cbhcompoundbatch ALTER COLUMN id SET DEFAULT nextval('cbh_chembl_model_extension_cbhcompoundbatch_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY corsheaders_corsmodel ALTER COLUMN id SET DEFAULT nextval('corsheaders_corsmodel_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY django_content_type ALTER COLUMN id SET DEFAULT nextval('django_content_type_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY django_site ALTER COLUMN id SET DEFAULT nextval('django_site_id_seq'::regclass);


--
-- Name: molregno; Type: DEFAULT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY molecule_dictionary ALTER COLUMN molregno SET DEFAULT nextval('molregno_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY south_migrationhistory ALTER COLUMN id SET DEFAULT nextval('south_migrationhistory_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY tastypie_apiaccess ALTER COLUMN id SET DEFAULT nextval('tastypie_apiaccess_id_seq'::regclass);


--
-- Data for Name: action_type; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY action_type (action_type, description, parent_type) FROM stdin;
\.


--
-- Data for Name: activities; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY activities (activity_id, assay_id, doc_id, record_id, molregno, standard_relation, published_value, published_units, standard_value, standard_units, standard_flag, standard_type, activity_comment, published_type, data_validity_comment, potential_duplicate, published_relation, pchembl_value, bao_endpoint, uo_units, qudt_units) FROM stdin;
\.


--
-- Data for Name: activity_stds_lookup; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY activity_stds_lookup (std_act_id, standard_type, definition, standard_units, normal_range_min, normal_range_max) FROM stdin;
\.


--
-- Data for Name: assay_parameters; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY assay_parameters (assay_param_id, assay_id, parameter_type, parameter_value) FROM stdin;
\.


--
-- Data for Name: assay_type; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY assay_type (assay_type, assay_desc) FROM stdin;
\.


--
-- Data for Name: assays; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY assays (assay_id, doc_id, description, assay_type, assay_test_type, assay_category, assay_organism, assay_tax_id, assay_strain, assay_tissue, assay_cell_type, assay_subcellular_fraction, tid, relationship_type, confidence_score, curated_by, src_id, src_assay_id, chembl_id, cell_id, bao_format) FROM stdin;
\.


--
-- Data for Name: atc_classification; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY atc_classification (who_name, level1, level2, level3, level4, level5, who_id, level1_description, level2_description, level3_description, level4_description) FROM stdin;
\.


--
-- Data for Name: auth_group; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY auth_group (id, name) FROM stdin;
\.


--
-- Name: auth_group_id_seq; Type: SEQUENCE SET; Schema: public; Owner: chembl
--

SELECT pg_catalog.setval('auth_group_id_seq', 1, false);


--
-- Data for Name: auth_group_permissions; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY auth_group_permissions (id, group_id, permission_id) FROM stdin;
\.


--
-- Name: auth_group_permissions_id_seq; Type: SEQUENCE SET; Schema: public; Owner: chembl
--

SELECT pg_catalog.setval('auth_group_permissions_id_seq', 1, false);


--
-- Data for Name: auth_permission; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY auth_permission (id, name, content_type_id, codename) FROM stdin;
1	Can add content type	1	add_contenttype
2	Can change content type	1	change_contenttype
3	Can delete content type	1	delete_contenttype
4	Can add site	2	add_site
5	Can change site	2	change_site
6	Can delete site	2	delete_site
7	Can add permission	5	add_permission
8	Can change permission	5	change_permission
9	Can delete permission	5	delete_permission
10	Can add group	3	add_group
11	Can change group	3	change_group
12	Can delete group	3	delete_group
13	Can add user	4	add_user
14	Can change user	4	change_user
15	Can delete user	4	delete_user
16	Can add cors model	6	add_corsmodel
17	Can change cors model	6	change_corsmodel
18	Can delete cors model	6	delete_corsmodel
19	Can add migration history	7	add_migrationhistory
20	Can change migration history	7	change_migrationhistory
21	Can delete migration history	7	delete_migrationhistory
22	Can add chembl id lookup	67	add_chemblidlookup
23	Can change chembl id lookup	67	change_chemblidlookup
24	Can delete chembl id lookup	67	delete_chemblidlookup
25	Can add version	68	add_version
26	Can change version	68	change_version
27	Can delete version	68	delete_version
28	Can add inchi errors	66	add_inchierrors
29	Can change inchi errors	66	change_inchierrors
30	Can delete inchi errors	66	delete_inchierrors
31	Can add image errors	8	add_imageerrors
32	Can change image errors	8	change_imageerrors
33	Can delete image errors	8	delete_imageerrors
34	Can add django cheat sheet	46	add_djangocheatsheet
35	Can change django cheat sheet	46	change_djangocheatsheet
36	Can delete django cheat sheet	46	delete_djangocheatsheet
37	Can add sdf	42	add_sdf
38	Can change sdf	42	change_sdf
39	Can delete sdf	42	delete_sdf
40	Can add journals	69	add_journals
41	Can change journals	69	change_journals
42	Can delete journals	69	delete_journals
43	Can add journal articles	70	add_journalarticles
44	Can change journal articles	70	change_journalarticles
45	Can delete journal articles	70	delete_journalarticles
46	Can add docs	71	add_docs
47	Can change docs	71	change_docs
48	Can delete docs	71	delete_docs
49	Can add source	72	add_source
50	Can change source	72	change_source
51	Can delete source	72	delete_source
52	Can add molecule dictionary	73	add_moleculedictionary
53	Can change molecule dictionary	73	change_moleculedictionary
54	Can delete molecule dictionary	73	delete_moleculedictionary
55	Can add compound records	74	add_compoundrecords
56	Can change compound records	74	change_compoundrecords
57	Can delete compound records	74	delete_compoundrecords
58	Can add compound properties	75	add_compoundproperties
59	Can change compound properties	75	change_compoundproperties
60	Can delete compound properties	75	delete_compoundproperties
61	Can add molecule hierarchy	76	add_moleculehierarchy
62	Can change molecule hierarchy	76	change_moleculehierarchy
63	Can delete molecule hierarchy	76	delete_moleculehierarchy
64	Can add research stem	77	add_researchstem
65	Can change research stem	77	change_researchstem
66	Can delete research stem	77	delete_researchstem
67	Can add research companies	78	add_researchcompanies
68	Can change research companies	78	change_researchcompanies
69	Can delete research companies	78	delete_researchcompanies
70	Can add molecule synonyms	79	add_moleculesynonyms
71	Can change molecule synonyms	79	change_moleculesynonyms
72	Can delete molecule synonyms	79	delete_moleculesynonyms
73	Can add biotherapeutics	80	add_biotherapeutics
74	Can change biotherapeutics	80	change_biotherapeutics
75	Can delete biotherapeutics	80	delete_biotherapeutics
76	Can add bio component sequences	81	add_biocomponentsequences
77	Can change bio component sequences	81	change_biocomponentsequences
78	Can delete bio component sequences	81	delete_biocomponentsequences
79	Can add biotherapeutic components	82	add_biotherapeuticcomponents
80	Can change biotherapeutic components	82	change_biotherapeuticcomponents
81	Can delete biotherapeutic components	82	delete_biotherapeuticcomponents
82	Can add compound images	83	add_compoundimages
83	Can change compound images	83	change_compoundimages
84	Can delete compound images	83	delete_compoundimages
85	Can add compound mols	84	add_compoundmols
86	Can change compound mols	84	change_compoundmols
87	Can delete compound mols	84	delete_compoundmols
88	Can add compound structures	85	add_compoundstructures
89	Can change compound structures	85	change_compoundstructures
90	Can delete compound structures	85	delete_compoundstructures
91	Can add target type	86	add_targettype
92	Can change target type	86	change_targettype
93	Can delete target type	86	delete_targettype
94	Can add target dictionary	87	add_targetdictionary
95	Can change target dictionary	87	change_targetdictionary
96	Can delete target dictionary	87	delete_targetdictionary
97	Can add component sequences	88	add_componentsequences
98	Can change component sequences	88	change_componentsequences
99	Can delete component sequences	88	delete_componentsequences
100	Can add component synonyms	89	add_componentsynonyms
101	Can change component synonyms	89	change_componentsynonyms
102	Can delete component synonyms	89	delete_componentsynonyms
103	Can add target components	90	add_targetcomponents
104	Can change target components	90	change_targetcomponents
105	Can delete target components	90	delete_targetcomponents
106	Can add organism class	91	add_organismclass
107	Can change organism class	91	change_organismclass
108	Can delete organism class	91	delete_organismclass
109	Can add protein family classification	92	add_proteinfamilyclassification
110	Can change protein family classification	92	change_proteinfamilyclassification
111	Can delete protein family classification	92	delete_proteinfamilyclassification
112	Can add component class	93	add_componentclass
113	Can change component class	93	change_componentclass
114	Can delete component class	93	delete_componentclass
115	Can add target relations	94	add_targetrelations
116	Can change target relations	94	change_targetrelations
117	Can delete target relations	94	delete_targetrelations
118	Can add protein classification	95	add_proteinclassification
119	Can change protein classification	95	change_proteinclassification
120	Can delete protein classification	95	delete_proteinclassification
121	Can add protein class synonyms	96	add_proteinclasssynonyms
122	Can change protein class synonyms	96	change_proteinclasssynonyms
123	Can delete protein class synonyms	96	delete_proteinclasssynonyms
124	Can add binding sites	97	add_bindingsites
125	Can change binding sites	97	change_bindingsites
126	Can delete binding sites	97	delete_bindingsites
127	Can add domains	98	add_domains
128	Can change domains	98	change_domains
129	Can delete domains	98	delete_domains
130	Can add component domains	99	add_componentdomains
131	Can change component domains	99	change_componentdomains
132	Can delete component domains	99	delete_componentdomains
133	Can add site components	100	add_sitecomponents
134	Can change site components	100	change_sitecomponents
135	Can delete site components	100	delete_sitecomponents
136	Can add relationship type	101	add_relationshiptype
137	Can change relationship type	101	change_relationshiptype
138	Can delete relationship type	101	delete_relationshiptype
139	Can add confidence score lookup	102	add_confidencescorelookup
140	Can change confidence score lookup	102	change_confidencescorelookup
141	Can delete confidence score lookup	102	delete_confidencescorelookup
142	Can add curation lookup	103	add_curationlookup
143	Can change curation lookup	103	change_curationlookup
144	Can delete curation lookup	103	delete_curationlookup
145	Can add assay type	104	add_assaytype
146	Can change assay type	104	change_assaytype
147	Can delete assay type	104	delete_assaytype
148	Can add assays	105	add_assays
149	Can change assays	105	change_assays
150	Can delete assays	105	delete_assays
151	Can add data validity lookup	106	add_datavaliditylookup
152	Can change data validity lookup	106	change_datavaliditylookup
153	Can delete data validity lookup	106	delete_datavaliditylookup
154	Can add parameter type	107	add_parametertype
155	Can change parameter type	107	change_parametertype
156	Can delete parameter type	107	delete_parametertype
157	Can add assay parameters	108	add_assayparameters
158	Can change assay parameters	108	change_assayparameters
159	Can delete assay parameters	108	delete_assayparameters
160	Can add activities	109	add_activities
161	Can change activities	109	change_activities
162	Can delete activities	109	delete_activities
163	Can add activity stds lookup	110	add_activitystdslookup
164	Can change activity stds lookup	110	change_activitystdslookup
165	Can delete activity stds lookup	110	delete_activitystdslookup
166	Can add products	111	add_products
167	Can change products	111	change_products
168	Can delete products	111	delete_products
169	Can add formulations	112	add_formulations
170	Can change formulations	112	change_formulations
171	Can delete formulations	112	delete_formulations
172	Can add atc classification	113	add_atcclassification
173	Can change atc classification	113	change_atcclassification
174	Can delete atc classification	113	delete_atcclassification
175	Can add defined daily dose	114	add_defineddailydose
176	Can change defined daily dose	114	change_defineddailydose
177	Can delete defined daily dose	114	delete_defineddailydose
178	Can add usan stems	115	add_usanstems
179	Can change usan stems	115	change_usanstems
180	Can delete usan stems	115	delete_usanstems
181	Can add molecule atc classification	116	add_moleculeatcclassification
182	Can change molecule atc classification	116	change_moleculeatcclassification
183	Can delete molecule atc classification	116	delete_moleculeatcclassification
184	Can add predicted binding domains	117	add_predictedbindingdomains
185	Can change predicted binding domains	117	change_predictedbindingdomains
186	Can delete predicted binding domains	117	delete_predictedbindingdomains
187	Can add ligand eff	118	add_ligandeff
188	Can change ligand eff	118	change_ligandeff
189	Can delete ligand eff	118	delete_ligandeff
190	Can add action type	119	add_actiontype
191	Can change action type	119	change_actiontype
192	Can delete action type	119	delete_actiontype
193	Can add drug mechanism	120	add_drugmechanism
194	Can change drug mechanism	120	change_drugmechanism
195	Can delete drug mechanism	120	delete_drugmechanism
196	Can add mechanism refs	121	add_mechanismrefs
197	Can change mechanism refs	121	change_mechanismrefs
198	Can delete mechanism refs	121	delete_mechanismrefs
199	Can add session	122	add_session
200	Can change session	122	change_session
201	Can delete session	122	delete_session
202	Can add cbh compound batch	123	add_cbhcompoundbatch
203	Can change cbh compound batch	123	change_cbhcompoundbatch
204	Can delete cbh compound batch	123	delete_cbhcompoundbatch
\.


--
-- Name: auth_permission_id_seq; Type: SEQUENCE SET; Schema: public; Owner: chembl
--

SELECT pg_catalog.setval('auth_permission_id_seq', 204, true);


--
-- Data for Name: auth_user; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY auth_user (id, password, last_login, is_superuser, username, first_name, last_name, email, is_staff, is_active, date_joined) FROM stdin;
1	pbkdf2_sha256$10000$ud7qhRa8UW8I$B/DmurvmEOneL6A53jitdHJ1n2b8rzPQlC0LbzyaGOY=	2014-12-17 11:40:49.718393+00	t	andy			strets123@gmail.com	t	t	2014-12-12 12:57:08.263386+00
\.


--
-- Data for Name: auth_user_groups; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY auth_user_groups (id, user_id, group_id) FROM stdin;
\.


--
-- Name: auth_user_groups_id_seq; Type: SEQUENCE SET; Schema: public; Owner: chembl
--

SELECT pg_catalog.setval('auth_user_groups_id_seq', 1, false);


--
-- Name: auth_user_id_seq; Type: SEQUENCE SET; Schema: public; Owner: chembl
--

SELECT pg_catalog.setval('auth_user_id_seq', 1, true);


--
-- Data for Name: auth_user_user_permissions; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY auth_user_user_permissions (id, user_id, permission_id) FROM stdin;
\.


--
-- Name: auth_user_user_permissions_id_seq; Type: SEQUENCE SET; Schema: public; Owner: chembl
--

SELECT pg_catalog.setval('auth_user_user_permissions_id_seq', 1, false);


--
-- Data for Name: binding_sites; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY binding_sites (site_id, site_name, tid) FROM stdin;
\.


--
-- Data for Name: bio_component_sequences; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY bio_component_sequences (component_id, component_type, description, sequence, sequence_md5sum, tax_id, organism) FROM stdin;
\.


--
-- Data for Name: biotherapeutic_components; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY biotherapeutic_components (biocomp_id, molregno, component_id) FROM stdin;
\.


--
-- Data for Name: biotherapeutics; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY biotherapeutics (molregno, description) FROM stdin;
\.


--
-- Data for Name: cbh_chembl_model_extension_cbhcompoundbatch; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY cbh_chembl_model_extension_cbhcompoundbatch (id, ctab, editable_by, viewable_by, related_molregno_id, warnings, custom_fields) FROM stdin;
\.


--
-- Name: cbh_chembl_model_extension_cbhcompoundbatch_id_seq; Type: SEQUENCE SET; Schema: public; Owner: chembl
--

SELECT pg_catalog.setval('cbh_chembl_model_extension_cbhcompoundbatch_id_seq', 1, false);


--
-- Data for Name: cell_dictionary; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY cell_dictionary (cell_id, cell_name, cell_description, cell_source_tissue, cell_source_organism, cell_source_tax_id, clo_id, efo_id, cellosaurus_id) FROM stdin;
\.


--
-- Data for Name: chembl_id_lookup; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY chembl_id_lookup (chembl_id, entity_type, entity_id, status) FROM stdin;
sdfkfsdfjsdf	COMPOUND	1	ACTIVE
\.


--
-- Data for Name: component_class; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY component_class (component_id, protein_class_id, comp_class_id) FROM stdin;
\.


--
-- Data for Name: component_domains; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY component_domains (compd_id, domain_id, component_id, start_position, end_position) FROM stdin;
\.


--
-- Data for Name: component_sequences; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY component_sequences (component_id, component_type, accession, sequence, sequence_md5sum, description, tax_id, organism, db_source, db_version) FROM stdin;
\.


--
-- Data for Name: component_synonyms; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY component_synonyms (compsyn_id, component_id, component_synonym, syn_type) FROM stdin;
\.


--
-- Data for Name: compound_mols; Type: TABLE DATA; Schema: public; Owner: vagrant
--

COPY compound_mols (molregno, ctab) FROM stdin;
\.


--
-- Data for Name: compound_properties; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY compound_properties (molregno, mw_freebase, alogp, hba, hbd, psa, rtb, ro3_pass, num_ro5_violations, med_chem_friendly, acd_most_apka, acd_most_bpka, acd_logp, acd_logd, molecular_species, full_mwt, aromatic_rings, heavy_atoms, num_alerts, qed_weighted, mw_monoisotopic, full_molformula) FROM stdin;
\.


--
-- Data for Name: compound_records; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY compound_records (record_id, molregno, doc_id, compound_key, compound_name, src_id, src_compound_id) FROM stdin;
\.


--
-- Data for Name: compound_structures; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY compound_structures (molregno, molfile, standard_inchi, standard_inchi_key, canonical_smiles) FROM stdin;
\.


--
-- Data for Name: confidence_score_lookup; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY confidence_score_lookup (confidence_score, description, target_mapping) FROM stdin;
\.


--
-- Data for Name: corsheaders_corsmodel; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY corsheaders_corsmodel (id, cors) FROM stdin;
\.


--
-- Name: corsheaders_corsmodel_id_seq; Type: SEQUENCE SET; Schema: public; Owner: chembl
--

SELECT pg_catalog.setval('corsheaders_corsmodel_id_seq', 1, false);


--
-- Data for Name: curation_lookup; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY curation_lookup (curated_by, description) FROM stdin;
\.


--
-- Data for Name: data_validity_lookup; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY data_validity_lookup (data_validity_comment, description) FROM stdin;
\.


--
-- Data for Name: defined_daily_dose; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY defined_daily_dose (atc_code, ddd_value, ddd_units, ddd_admr, ddd_comment, ddd_id) FROM stdin;
\.


--
-- Data for Name: django_content_type; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY django_content_type (id, name, app_label, model) FROM stdin;
1	content type	contenttypes	contenttype
2	site	sites	site
3	group	auth	group
4	user	auth	user
5	permission	auth	permission
6	cors model	corsheaders	corsmodel
7	migration history	south	migrationhistory
8	image errors	chembl_business_model	imageerrors
9	component synonyms	chembl_business_model	componentsynonyms
10	target dictionary	chembl_business_model	targetdictionary
11	component sequences	chembl_business_model	componentsequences
12	compound records	chembl_business_model	compoundrecords
13	action type	chembl_business_model	actiontype
14	target components	chembl_business_model	targetcomponents
15	confidence score lookup	chembl_business_model	confidencescorelookup
16	molecule atc classification	chembl_business_model	moleculeatcclassification
17	component domains	chembl_business_model	componentdomains
18	component class	chembl_business_model	componentclass
19	assay parameters	chembl_business_model	assayparameters
20	research companies	chembl_business_model	researchcompanies
21	journal articles	chembl_business_model	journalarticles
22	biotherapeutics	chembl_business_model	biotherapeutics
23	site components	chembl_business_model	sitecomponents
24	biotherapeutic components	chembl_business_model	biotherapeuticcomponents
25	curation lookup	chembl_business_model	curationlookup
26	protein classification	chembl_business_model	proteinclassification
27	atc classification	chembl_business_model	atcclassification
28	formulations	chembl_business_model	formulations
29	activities	chembl_business_model	activities
30	source	chembl_business_model	source
31	version	chembl_business_model	version
32	activity stds lookup	chembl_business_model	activitystdslookup
33	products	chembl_business_model	products
34	compound structures	chembl_business_model	compoundstructures
35	predicted binding domains	chembl_business_model	predictedbindingdomains
36	molecule synonyms	chembl_business_model	moleculesynonyms
37	assay type	chembl_business_model	assaytype
38	target type	chembl_business_model	targettype
39	docs	chembl_business_model	docs
40	compound mols	chembl_business_model	compoundmols
41	bio component sequences	chembl_business_model	biocomponentsequences
42	sdf	chembl_business_model	sdf
43	molecule hierarchy	chembl_business_model	moleculehierarchy
44	protein family classification	chembl_business_model	proteinfamilyclassification
45	drug mechanism	chembl_business_model	drugmechanism
46	django cheat sheet	chembl_business_model	djangocheatsheet
47	usan stems	chembl_business_model	usanstems
48	relationship type	chembl_business_model	relationshiptype
49	journals	chembl_business_model	journals
50	defined daily dose	chembl_business_model	defineddailydose
51	chembl id lookup	chembl_business_model	chemblidlookup
52	compound images	chembl_business_model	compoundimages
53	parameter type	chembl_business_model	parametertype
54	target relations	chembl_business_model	targetrelations
55	research stem	chembl_business_model	researchstem
56	binding sites	chembl_business_model	bindingsites
57	data validity lookup	chembl_business_model	datavaliditylookup
58	ligand eff	chembl_business_model	ligandeff
59	compound properties	chembl_business_model	compoundproperties
60	mechanism refs	chembl_business_model	mechanismrefs
61	organism class	chembl_business_model	organismclass
62	molecule dictionary	chembl_business_model	moleculedictionary
63	protein class synonyms	chembl_business_model	proteinclasssynonyms
64	assays	chembl_business_model	assays
65	domains	chembl_business_model	domains
66	inchi errors	chembl_business_model	inchierrors
67	chembl id lookup	chembl_core_model	chemblidlookup
68	version	chembl_core_model	version
69	journals	chembl_core_model	journals
70	journal articles	chembl_core_model	journalarticles
71	docs	chembl_core_model	docs
72	source	chembl_core_model	source
73	molecule dictionary	chembl_core_model	moleculedictionary
74	compound records	chembl_core_model	compoundrecords
75	compound properties	chembl_core_model	compoundproperties
76	molecule hierarchy	chembl_core_model	moleculehierarchy
77	research stem	chembl_core_model	researchstem
78	research companies	chembl_core_model	researchcompanies
79	molecule synonyms	chembl_core_model	moleculesynonyms
80	biotherapeutics	chembl_core_model	biotherapeutics
81	bio component sequences	chembl_core_model	biocomponentsequences
82	biotherapeutic components	chembl_core_model	biotherapeuticcomponents
83	compound images	chembl_core_model	compoundimages
84	compound mols	chembl_core_model	compoundmols
85	compound structures	chembl_core_model	compoundstructures
86	target type	chembl_core_model	targettype
87	target dictionary	chembl_core_model	targetdictionary
88	component sequences	chembl_core_model	componentsequences
89	component synonyms	chembl_core_model	componentsynonyms
90	target components	chembl_core_model	targetcomponents
91	organism class	chembl_core_model	organismclass
92	protein family classification	chembl_core_model	proteinfamilyclassification
93	component class	chembl_core_model	componentclass
94	target relations	chembl_core_model	targetrelations
95	protein classification	chembl_core_model	proteinclassification
96	protein class synonyms	chembl_core_model	proteinclasssynonyms
97	binding sites	chembl_core_model	bindingsites
98	domains	chembl_core_model	domains
99	component domains	chembl_core_model	componentdomains
100	site components	chembl_core_model	sitecomponents
101	relationship type	chembl_core_model	relationshiptype
102	confidence score lookup	chembl_core_model	confidencescorelookup
103	curation lookup	chembl_core_model	curationlookup
104	assay type	chembl_core_model	assaytype
105	assays	chembl_core_model	assays
106	data validity lookup	chembl_core_model	datavaliditylookup
107	parameter type	chembl_core_model	parametertype
108	assay parameters	chembl_core_model	assayparameters
109	activities	chembl_core_model	activities
110	activity stds lookup	chembl_core_model	activitystdslookup
111	products	chembl_core_model	products
112	formulations	chembl_core_model	formulations
113	atc classification	chembl_core_model	atcclassification
114	defined daily dose	chembl_core_model	defineddailydose
115	usan stems	chembl_core_model	usanstems
116	molecule atc classification	chembl_core_model	moleculeatcclassification
117	predicted binding domains	chembl_core_model	predictedbindingdomains
118	ligand eff	chembl_core_model	ligandeff
119	action type	chembl_core_model	actiontype
120	drug mechanism	chembl_core_model	drugmechanism
121	mechanism refs	chembl_core_model	mechanismrefs
122	session	sessions	session
123	cbh compound batch	cbh_chembl_model_extension	cbhcompoundbatch
\.


--
-- Name: django_content_type_id_seq; Type: SEQUENCE SET; Schema: public; Owner: chembl
--

SELECT pg_catalog.setval('django_content_type_id_seq', 123, true);


--
-- Data for Name: django_session; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY django_session (session_key, session_data, expire_date) FROM stdin;
xd8ql55ksah6xnr7e8ov5nrgqocvj7s7	YTgxYTFlZjY0NmJkOTFlNmY3MDdhNjIxNmNhZmYwNTU5ZDU1NzBlZjqAAn1xAShVEl9hdXRoX3VzZXJfYmFja2VuZFUpZGphbmdvLmNvbnRyaWIuYXV0aC5iYWNrZW5kcy5Nb2RlbEJhY2tlbmRxAlUNX2F1dGhfdXNlcl9pZEsBdS4=	2014-12-30 11:42:57.994728+00
b374x2jzenwoxzdvqgpy5x3456ivhwpk	NTE0MWUyMjE4MWU5YjRmMDE3NTg0N2NlNDcwMzRjMGE1NWFiMTI0MTqAAn1xAShVEl9hdXRoX3VzZXJfYmFja2VuZHECVSlkamFuZ28uY29udHJpYi5hdXRoLmJhY2tlbmRzLk1vZGVsQmFja2VuZHEDVQ1fYXV0aF91c2VyX2lkcQRLAXUu	2014-12-31 11:40:49.72499+00
\.


--
-- Data for Name: django_site; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY django_site (id, domain, name) FROM stdin;
1	example.com	example.com
\.


--
-- Name: django_site_id_seq; Type: SEQUENCE SET; Schema: public; Owner: chembl
--

SELECT pg_catalog.setval('django_site_id_seq', 1, true);


--
-- Data for Name: docs; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY docs (doc_id, journal, year, volume, issue, first_page, last_page, pubmed_id, doi, chembl_id, title, doc_type, authors, abstract) FROM stdin;
\.


--
-- Data for Name: domains; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY domains (domain_id, domain_type, source_domain_id, domain_name, domain_description) FROM stdin;
\.


--
-- Data for Name: drug_mechanism; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY drug_mechanism (mec_id, record_id, molregno, mechanism_of_action, tid, site_id, action_type, direct_interaction, molecular_mechanism, disease_efficacy, mechanism_comment, selectivity_comment, binding_site_comment) FROM stdin;
\.


--
-- Data for Name: formulations; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY formulations (product_id, ingredient, strength, record_id, molregno, formulation_id) FROM stdin;
\.


--
-- Data for Name: fps2_rdkit; Type: TABLE DATA; Schema: public; Owner: vagrant
--

COPY fps2_rdkit (molregno, torsionbv, mfp2, ffp2, rdkfp, atombv, layeredfp, maccsfp) FROM stdin;
\.


--
-- Data for Name: fps_rdkit; Type: TABLE DATA; Schema: public; Owner: vagrant
--

COPY fps_rdkit (molregno, torsionbv, mfp2, ffp2) FROM stdin;
\.


--
-- Data for Name: ligand_eff; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY ligand_eff (activity_id, bei, sei, le, lle) FROM stdin;
\.


--
-- Data for Name: mechanism_refs; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY mechanism_refs (mecref_id, mec_id, ref_type, ref_id, ref_url) FROM stdin;
\.


--
-- Data for Name: molecule_atc_classification; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY molecule_atc_classification (mol_atc_id, level5, molregno) FROM stdin;
\.


--
-- Data for Name: molecule_dictionary; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY molecule_dictionary (molregno, pref_name, chembl_id, max_phase, therapeutic_flag, dosed_ingredient, structure_type, chebi_par_id, molecule_type, first_approval, oral, parenteral, topical, black_box_warning, natural_product, first_in_class, chirality, prodrug, inorganic_flag, usan_year, availability_type, usan_stem, polymer_flag, usan_substem, usan_stem_definition, indication_class, structure_key, nomerge_reason, chebi_id, insert_date, molfile_update, downgrade_reason, checked_by, replacement_mrn, exclude, downgraded, nomerge) FROM stdin;
6	\N	sdfkfsdfjsdf	0	0	0	MOL	\N	\N	\N	0	0	0	0	-1	-1	-1	-1	0	\N	\N	\N	\N	\N	\N	\N	\N	\N	\N	\N	\N	\N	\N	\N	\N	\N	\N
\.


--
-- Data for Name: molecule_hierarchy; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY molecule_hierarchy (molregno, parent_molregno, active_molregno) FROM stdin;
\.


--
-- Data for Name: molecule_synonyms; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY molecule_synonyms (molregno, syn_type, molsyn_id, res_stem_id, synonyms) FROM stdin;
\.


--
-- Name: molregno_id_seq; Type: SEQUENCE SET; Schema: public; Owner: chembl
--

SELECT pg_catalog.setval('molregno_id_seq', 6, true);


--
-- Data for Name: octmp_summary; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY octmp_summary (id, table_name, table_created, query_md5) FROM stdin;
\.


--
-- Name: octmp_summary_id_seq; Type: SEQUENCE SET; Schema: public; Owner: chembl
--

SELECT pg_catalog.setval('octmp_summary_id_seq', 1, false);


--
-- Data for Name: organism_class; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY organism_class (oc_id, tax_id, l1, l2, l3) FROM stdin;
\.


--
-- Data for Name: parameter_type; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY parameter_type (parameter_type, description) FROM stdin;
\.


--
-- Data for Name: predicted_binding_domains; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY predicted_binding_domains (predbind_id, activity_id, site_id, prediction_method, confidence) FROM stdin;
\.


--
-- Data for Name: products; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY products (dosage_form, route, trade_name, approval_date, ad_type, oral, topical, parenteral, black_box_warning, applicant_full_name, innovator_company, product_id, nda_type) FROM stdin;
\.


--
-- Data for Name: protein_class_synonyms; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY protein_class_synonyms (protclasssyn_id, protein_class_id, protein_class_synonym, syn_type) FROM stdin;
\.


--
-- Data for Name: protein_classification; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY protein_classification (protein_class_id, parent_id, pref_name, short_name, protein_class_desc, definition, class_level) FROM stdin;
\.


--
-- Data for Name: protein_family_classification; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY protein_family_classification (protein_class_id, protein_class_desc, l1, l2, l3, l4, l5, l6, l7, l8) FROM stdin;
\.


--
-- Data for Name: relationship_type; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY relationship_type (relationship_type, relationship_desc) FROM stdin;
\.


--
-- Data for Name: research_companies; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY research_companies (co_stem_id, res_stem_id, company, country, previous_company) FROM stdin;
\.


--
-- Data for Name: research_stem; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY research_stem (res_stem_id, research_stem) FROM stdin;
\.


--
-- Data for Name: site_components; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY site_components (sitecomp_id, site_id, component_id, domain_id, site_residues) FROM stdin;
\.


--
-- Data for Name: source; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY source (src_id, src_description, src_short_name) FROM stdin;
\.


--
-- Data for Name: south_migrationhistory; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY south_migrationhistory (id, app_name, migration, applied) FROM stdin;
12	cbh_chembl_model_extension	0001_initial	2014-12-16 14:13:09.796139+00
13	cbh_chembl_model_extension	0002_auto__chg_field_cbhcompoundbatch_related_molregno	2014-12-16 14:13:09.83201+00
\.


--
-- Name: south_migrationhistory_id_seq; Type: SEQUENCE SET; Schema: public; Owner: chembl
--

SELECT pg_catalog.setval('south_migrationhistory_id_seq', 13, true);


--
-- Data for Name: target_components; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY target_components (tid, component_id, targcomp_id, homologue) FROM stdin;
\.


--
-- Data for Name: target_dictionary; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY target_dictionary (tid, target_type, pref_name, tax_id, organism, chembl_id, species_group_flag) FROM stdin;
\.


--
-- Data for Name: target_relations; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY target_relations (tid, relationship, related_tid, targrel_id) FROM stdin;
\.


--
-- Data for Name: target_type; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY target_type (target_type, target_desc, parent_type) FROM stdin;
\.


--
-- Data for Name: tastypie_apiaccess; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY tastypie_apiaccess (id, identifier, url, request_method, accessed) FROM stdin;
\.


--
-- Name: tastypie_apiaccess_id_seq; Type: SEQUENCE SET; Schema: public; Owner: chembl
--

SELECT pg_catalog.setval('tastypie_apiaccess_id_seq', 1, false);


--
-- Data for Name: usan_stems; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY usan_stems (usan_stem_id, stem, subgroup, annotation, stem_class, major_class, who_extra) FROM stdin;
\.


--
-- Data for Name: version; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY version (name, creation_date, comments) FROM stdin;
\.


--
-- Data for Name: ws_cache; Type: TABLE DATA; Schema: public; Owner: chembl
--

COPY ws_cache (cache_key, value, expires) FROM stdin;
\.


--
-- Name: action_type_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY action_type
    ADD CONSTRAINT action_type_pkey PRIMARY KEY (action_type);


--
-- Name: activities_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY activities
    ADD CONSTRAINT activities_pkey PRIMARY KEY (activity_id);


--
-- Name: activity_stds_lookup_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY activity_stds_lookup
    ADD CONSTRAINT activity_stds_lookup_pkey PRIMARY KEY (std_act_id);


--
-- Name: activity_stds_lookup_standard_type_standard_units_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY activity_stds_lookup
    ADD CONSTRAINT activity_stds_lookup_standard_type_standard_units_key UNIQUE (standard_type, standard_units);


--
-- Name: assay_parameters_assay_id_parameter_type_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY assay_parameters
    ADD CONSTRAINT assay_parameters_assay_id_parameter_type_key UNIQUE (assay_id, parameter_type);


--
-- Name: assay_parameters_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY assay_parameters
    ADD CONSTRAINT assay_parameters_pkey PRIMARY KEY (assay_param_id);


--
-- Name: assay_type_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY assay_type
    ADD CONSTRAINT assay_type_pkey PRIMARY KEY (assay_type);


--
-- Name: assays_chembl_id_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY assays
    ADD CONSTRAINT assays_chembl_id_key UNIQUE (chembl_id);


--
-- Name: assays_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY assays
    ADD CONSTRAINT assays_pkey PRIMARY KEY (assay_id);


--
-- Name: atc_classification_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY atc_classification
    ADD CONSTRAINT atc_classification_pkey PRIMARY KEY (level5);


--
-- Name: auth_group_name_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY auth_group
    ADD CONSTRAINT auth_group_name_key UNIQUE (name);


--
-- Name: auth_group_permissions_group_id_permission_id_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY auth_group_permissions
    ADD CONSTRAINT auth_group_permissions_group_id_permission_id_key UNIQUE (group_id, permission_id);


--
-- Name: auth_group_permissions_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY auth_group_permissions
    ADD CONSTRAINT auth_group_permissions_pkey PRIMARY KEY (id);


--
-- Name: auth_group_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY auth_group
    ADD CONSTRAINT auth_group_pkey PRIMARY KEY (id);


--
-- Name: auth_permission_content_type_id_codename_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY auth_permission
    ADD CONSTRAINT auth_permission_content_type_id_codename_key UNIQUE (content_type_id, codename);


--
-- Name: auth_permission_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY auth_permission
    ADD CONSTRAINT auth_permission_pkey PRIMARY KEY (id);


--
-- Name: auth_user_groups_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY auth_user_groups
    ADD CONSTRAINT auth_user_groups_pkey PRIMARY KEY (id);


--
-- Name: auth_user_groups_user_id_group_id_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY auth_user_groups
    ADD CONSTRAINT auth_user_groups_user_id_group_id_key UNIQUE (user_id, group_id);


--
-- Name: auth_user_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY auth_user
    ADD CONSTRAINT auth_user_pkey PRIMARY KEY (id);


--
-- Name: auth_user_user_permissions_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY auth_user_user_permissions
    ADD CONSTRAINT auth_user_user_permissions_pkey PRIMARY KEY (id);


--
-- Name: auth_user_user_permissions_user_id_permission_id_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY auth_user_user_permissions
    ADD CONSTRAINT auth_user_user_permissions_user_id_permission_id_key UNIQUE (user_id, permission_id);


--
-- Name: auth_user_username_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY auth_user
    ADD CONSTRAINT auth_user_username_key UNIQUE (username);


--
-- Name: binding_sites_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY binding_sites
    ADD CONSTRAINT binding_sites_pkey PRIMARY KEY (site_id);


--
-- Name: bio_component_sequences_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY bio_component_sequences
    ADD CONSTRAINT bio_component_sequences_pkey PRIMARY KEY (component_id);


--
-- Name: biotherapeutic_components_molregno_component_id_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY biotherapeutic_components
    ADD CONSTRAINT biotherapeutic_components_molregno_component_id_key UNIQUE (molregno, component_id);


--
-- Name: biotherapeutic_components_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY biotherapeutic_components
    ADD CONSTRAINT biotherapeutic_components_pkey PRIMARY KEY (biocomp_id);


--
-- Name: biotherapeutics_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY biotherapeutics
    ADD CONSTRAINT biotherapeutics_pkey PRIMARY KEY (molregno);


--
-- Name: cbh_chembl_model_extension_cbhcompoundbatch_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY cbh_chembl_model_extension_cbhcompoundbatch
    ADD CONSTRAINT cbh_chembl_model_extension_cbhcompoundbatch_pkey PRIMARY KEY (id);


--
-- Name: cell_dictionary_cell_name_cell_source_tax_id_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY cell_dictionary
    ADD CONSTRAINT cell_dictionary_cell_name_cell_source_tax_id_key UNIQUE (cell_name, cell_source_tax_id);


--
-- Name: cell_dictionary_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY cell_dictionary
    ADD CONSTRAINT cell_dictionary_pkey PRIMARY KEY (cell_id);


--
-- Name: chembl_id_lookup_entity_id_entity_type_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY chembl_id_lookup
    ADD CONSTRAINT chembl_id_lookup_entity_id_entity_type_key UNIQUE (entity_id, entity_type);


--
-- Name: chembl_id_lookup_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY chembl_id_lookup
    ADD CONSTRAINT chembl_id_lookup_pkey PRIMARY KEY (chembl_id);


--
-- Name: component_class_component_id_protein_class_id_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY component_class
    ADD CONSTRAINT component_class_component_id_protein_class_id_key UNIQUE (component_id, protein_class_id);


--
-- Name: component_class_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY component_class
    ADD CONSTRAINT component_class_pkey PRIMARY KEY (comp_class_id);


--
-- Name: component_domains_domain_id_component_id_start_position_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY component_domains
    ADD CONSTRAINT component_domains_domain_id_component_id_start_position_key UNIQUE (domain_id, component_id, start_position);


--
-- Name: component_domains_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY component_domains
    ADD CONSTRAINT component_domains_pkey PRIMARY KEY (compd_id);


--
-- Name: component_sequences_accession_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY component_sequences
    ADD CONSTRAINT component_sequences_accession_key UNIQUE (accession);


--
-- Name: component_sequences_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY component_sequences
    ADD CONSTRAINT component_sequences_pkey PRIMARY KEY (component_id);


--
-- Name: component_synonyms_component_id_component_synonym_syn_type_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY component_synonyms
    ADD CONSTRAINT component_synonyms_component_id_component_synonym_syn_type_key UNIQUE (component_id, component_synonym, syn_type);


--
-- Name: component_synonyms_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY component_synonyms
    ADD CONSTRAINT component_synonyms_pkey PRIMARY KEY (compsyn_id);


--
-- Name: compound_mols_pkey; Type: CONSTRAINT; Schema: public; Owner: vagrant; Tablespace: 
--

ALTER TABLE ONLY compound_mols
    ADD CONSTRAINT compound_mols_pkey PRIMARY KEY (molregno);


--
-- Name: compound_properties_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY compound_properties
    ADD CONSTRAINT compound_properties_pkey PRIMARY KEY (molregno);


--
-- Name: compound_records_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY compound_records
    ADD CONSTRAINT compound_records_pkey PRIMARY KEY (record_id);


--
-- Name: compound_structures_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY compound_structures
    ADD CONSTRAINT compound_structures_pkey PRIMARY KEY (molregno);


--
-- Name: confidence_score_lookup_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY confidence_score_lookup
    ADD CONSTRAINT confidence_score_lookup_pkey PRIMARY KEY (confidence_score);


--
-- Name: corsheaders_corsmodel_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY corsheaders_corsmodel
    ADD CONSTRAINT corsheaders_corsmodel_pkey PRIMARY KEY (id);


--
-- Name: curation_lookup_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY curation_lookup
    ADD CONSTRAINT curation_lookup_pkey PRIMARY KEY (curated_by);


--
-- Name: data_validity_lookup_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY data_validity_lookup
    ADD CONSTRAINT data_validity_lookup_pkey PRIMARY KEY (data_validity_comment);


--
-- Name: defined_daily_dose_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY defined_daily_dose
    ADD CONSTRAINT defined_daily_dose_pkey PRIMARY KEY (ddd_id);


--
-- Name: django_content_type_app_label_model_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY django_content_type
    ADD CONSTRAINT django_content_type_app_label_model_key UNIQUE (app_label, model);


--
-- Name: django_content_type_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY django_content_type
    ADD CONSTRAINT django_content_type_pkey PRIMARY KEY (id);


--
-- Name: django_session_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY django_session
    ADD CONSTRAINT django_session_pkey PRIMARY KEY (session_key);


--
-- Name: django_site_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY django_site
    ADD CONSTRAINT django_site_pkey PRIMARY KEY (id);


--
-- Name: docs_chembl_id_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY docs
    ADD CONSTRAINT docs_chembl_id_key UNIQUE (chembl_id);


--
-- Name: docs_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY docs
    ADD CONSTRAINT docs_pkey PRIMARY KEY (doc_id);


--
-- Name: docs_pubmed_id_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY docs
    ADD CONSTRAINT docs_pubmed_id_key UNIQUE (pubmed_id);


--
-- Name: domains_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY domains
    ADD CONSTRAINT domains_pkey PRIMARY KEY (domain_id);


--
-- Name: drug_mechanism_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY drug_mechanism
    ADD CONSTRAINT drug_mechanism_pkey PRIMARY KEY (mec_id);


--
-- Name: formulations_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY formulations
    ADD CONSTRAINT formulations_pkey PRIMARY KEY (formulation_id);


--
-- Name: formulations_record_id_product_id_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY formulations
    ADD CONSTRAINT formulations_record_id_product_id_key UNIQUE (record_id, product_id);


--
-- Name: fps2_rdkit_pkey; Type: CONSTRAINT; Schema: public; Owner: vagrant; Tablespace: 
--

ALTER TABLE ONLY fps2_rdkit
    ADD CONSTRAINT fps2_rdkit_pkey PRIMARY KEY (molregno);


--
-- Name: fps_rdkit_pkey; Type: CONSTRAINT; Schema: public; Owner: vagrant; Tablespace: 
--

ALTER TABLE ONLY fps_rdkit
    ADD CONSTRAINT fps_rdkit_pkey PRIMARY KEY (molregno);


--
-- Name: ligand_eff_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY ligand_eff
    ADD CONSTRAINT ligand_eff_pkey PRIMARY KEY (activity_id);


--
-- Name: mechanism_refs_mec_id_ref_type_ref_id_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY mechanism_refs
    ADD CONSTRAINT mechanism_refs_mec_id_ref_type_ref_id_key UNIQUE (mec_id, ref_type, ref_id);


--
-- Name: mechanism_refs_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY mechanism_refs
    ADD CONSTRAINT mechanism_refs_pkey PRIMARY KEY (mecref_id);


--
-- Name: molecule_atc_classification_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY molecule_atc_classification
    ADD CONSTRAINT molecule_atc_classification_pkey PRIMARY KEY (mol_atc_id);


--
-- Name: molecule_dictionary_chembl_id_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY molecule_dictionary
    ADD CONSTRAINT molecule_dictionary_chembl_id_key UNIQUE (chembl_id);


--
-- Name: molecule_dictionary_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY molecule_dictionary
    ADD CONSTRAINT molecule_dictionary_pkey PRIMARY KEY (molregno);


--
-- Name: molecule_hierarchy_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY molecule_hierarchy
    ADD CONSTRAINT molecule_hierarchy_pkey PRIMARY KEY (molregno);


--
-- Name: molecule_synonyms_molregno_synonyms_syn_type_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY molecule_synonyms
    ADD CONSTRAINT molecule_synonyms_molregno_synonyms_syn_type_key UNIQUE (molregno, synonyms, syn_type);


--
-- Name: molecule_synonyms_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY molecule_synonyms
    ADD CONSTRAINT molecule_synonyms_pkey PRIMARY KEY (molsyn_id);


--
-- Name: organism_class_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY organism_class
    ADD CONSTRAINT organism_class_pkey PRIMARY KEY (oc_id);


--
-- Name: organism_class_tax_id_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY organism_class
    ADD CONSTRAINT organism_class_tax_id_key UNIQUE (tax_id);


--
-- Name: parameter_type_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY parameter_type
    ADD CONSTRAINT parameter_type_pkey PRIMARY KEY (parameter_type);


--
-- Name: predicted_binding_domains_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY predicted_binding_domains
    ADD CONSTRAINT predicted_binding_domains_pkey PRIMARY KEY (predbind_id);


--
-- Name: products_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY products
    ADD CONSTRAINT products_pkey PRIMARY KEY (product_id);


--
-- Name: protein_class_synonyms_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY protein_class_synonyms
    ADD CONSTRAINT protein_class_synonyms_pkey PRIMARY KEY (protclasssyn_id);


--
-- Name: protein_classification_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY protein_classification
    ADD CONSTRAINT protein_classification_pkey PRIMARY KEY (protein_class_id);


--
-- Name: protein_family_classification_l1_l2_l3_l4_l5_l6_l7_l8_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY protein_family_classification
    ADD CONSTRAINT protein_family_classification_l1_l2_l3_l4_l5_l6_l7_l8_key UNIQUE (l1, l2, l3, l4, l5, l6, l7, l8);


--
-- Name: protein_family_classification_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY protein_family_classification
    ADD CONSTRAINT protein_family_classification_pkey PRIMARY KEY (protein_class_id);


--
-- Name: protein_family_classification_protein_class_desc_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY protein_family_classification
    ADD CONSTRAINT protein_family_classification_protein_class_desc_key UNIQUE (protein_class_desc);


--
-- Name: relationship_type_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY relationship_type
    ADD CONSTRAINT relationship_type_pkey PRIMARY KEY (relationship_type);


--
-- Name: research_companies_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY research_companies
    ADD CONSTRAINT research_companies_pkey PRIMARY KEY (co_stem_id);


--
-- Name: research_companies_res_stem_id_company_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY research_companies
    ADD CONSTRAINT research_companies_res_stem_id_company_key UNIQUE (res_stem_id, company);


--
-- Name: research_stem_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY research_stem
    ADD CONSTRAINT research_stem_pkey PRIMARY KEY (res_stem_id);


--
-- Name: research_stem_research_stem_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY research_stem
    ADD CONSTRAINT research_stem_research_stem_key UNIQUE (research_stem);


--
-- Name: site_components_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY site_components
    ADD CONSTRAINT site_components_pkey PRIMARY KEY (sitecomp_id);


--
-- Name: site_components_site_id_component_id_domain_id_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY site_components
    ADD CONSTRAINT site_components_site_id_component_id_domain_id_key UNIQUE (site_id, component_id, domain_id);


--
-- Name: source_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY source
    ADD CONSTRAINT source_pkey PRIMARY KEY (src_id);


--
-- Name: south_migrationhistory_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY south_migrationhistory
    ADD CONSTRAINT south_migrationhistory_pkey PRIMARY KEY (id);


--
-- Name: target_components_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY target_components
    ADD CONSTRAINT target_components_pkey PRIMARY KEY (targcomp_id);


--
-- Name: target_dictionary_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY target_dictionary
    ADD CONSTRAINT target_dictionary_pkey PRIMARY KEY (tid);


--
-- Name: target_relations_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY target_relations
    ADD CONSTRAINT target_relations_pkey PRIMARY KEY (targrel_id);


--
-- Name: target_type_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY target_type
    ADD CONSTRAINT target_type_pkey PRIMARY KEY (target_type);


--
-- Name: tastypie_apiaccess_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY tastypie_apiaccess
    ADD CONSTRAINT tastypie_apiaccess_pkey PRIMARY KEY (id);


--
-- Name: usan_stems_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY usan_stems
    ADD CONSTRAINT usan_stems_pkey PRIMARY KEY (usan_stem_id);


--
-- Name: usan_stems_stem_subgroup_key; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY usan_stems
    ADD CONSTRAINT usan_stems_stem_subgroup_key UNIQUE (stem, subgroup);


--
-- Name: version_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY version
    ADD CONSTRAINT version_pkey PRIMARY KEY (name);


--
-- Name: ws_cache_pkey; Type: CONSTRAINT; Schema: public; Owner: chembl; Tablespace: 
--

ALTER TABLE ONLY ws_cache
    ADD CONSTRAINT ws_cache_pkey PRIMARY KEY (cache_key);


--
-- Name: action_type_action_type_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX action_type_action_type_like ON action_type USING btree (action_type varchar_pattern_ops);


--
-- Name: activities_assay_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_assay_id ON activities USING btree (assay_id);


--
-- Name: activities_data_validity_comment; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_data_validity_comment ON activities USING btree (data_validity_comment);


--
-- Name: activities_data_validity_comment_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_data_validity_comment_like ON activities USING btree (data_validity_comment varchar_pattern_ops);


--
-- Name: activities_doc_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_doc_id ON activities USING btree (doc_id);


--
-- Name: activities_molregno; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_molregno ON activities USING btree (molregno);


--
-- Name: activities_pchembl_value; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_pchembl_value ON activities USING btree (pchembl_value);


--
-- Name: activities_published_relation; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_published_relation ON activities USING btree (published_relation);


--
-- Name: activities_published_relation_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_published_relation_like ON activities USING btree (published_relation varchar_pattern_ops);


--
-- Name: activities_published_type; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_published_type ON activities USING btree (published_type);


--
-- Name: activities_published_type_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_published_type_like ON activities USING btree (published_type varchar_pattern_ops);


--
-- Name: activities_published_units; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_published_units ON activities USING btree (published_units);


--
-- Name: activities_published_units_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_published_units_like ON activities USING btree (published_units varchar_pattern_ops);


--
-- Name: activities_published_value; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_published_value ON activities USING btree (published_value);


--
-- Name: activities_record_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_record_id ON activities USING btree (record_id);


--
-- Name: activities_standard_relation; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_standard_relation ON activities USING btree (standard_relation);


--
-- Name: activities_standard_relation_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_standard_relation_like ON activities USING btree (standard_relation varchar_pattern_ops);


--
-- Name: activities_standard_type; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_standard_type ON activities USING btree (standard_type);


--
-- Name: activities_standard_type_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_standard_type_like ON activities USING btree (standard_type varchar_pattern_ops);


--
-- Name: activities_standard_units; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_standard_units ON activities USING btree (standard_units);


--
-- Name: activities_standard_units_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_standard_units_like ON activities USING btree (standard_units varchar_pattern_ops);


--
-- Name: activities_standard_value; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX activities_standard_value ON activities USING btree (standard_value);


--
-- Name: assay_parameters_assay_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assay_parameters_assay_id ON assay_parameters USING btree (assay_id);


--
-- Name: assay_parameters_parameter_type; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assay_parameters_parameter_type ON assay_parameters USING btree (parameter_type);


--
-- Name: assay_parameters_parameter_type_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assay_parameters_parameter_type_like ON assay_parameters USING btree (parameter_type varchar_pattern_ops);


--
-- Name: assay_type_assay_type_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assay_type_assay_type_like ON assay_type USING btree (assay_type varchar_pattern_ops);


--
-- Name: assays_assay_type; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assays_assay_type ON assays USING btree (assay_type);


--
-- Name: assays_assay_type_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assays_assay_type_like ON assays USING btree (assay_type varchar_pattern_ops);


--
-- Name: assays_bao_format; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assays_bao_format ON assays USING btree (bao_format);


--
-- Name: assays_bao_format_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assays_bao_format_like ON assays USING btree (bao_format varchar_pattern_ops);


--
-- Name: assays_cell_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assays_cell_id ON assays USING btree (cell_id);


--
-- Name: assays_chembl_id_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assays_chembl_id_like ON assays USING btree (chembl_id varchar_pattern_ops);


--
-- Name: assays_confidence_score; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assays_confidence_score ON assays USING btree (confidence_score);


--
-- Name: assays_curated_by; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assays_curated_by ON assays USING btree (curated_by);


--
-- Name: assays_curated_by_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assays_curated_by_like ON assays USING btree (curated_by varchar_pattern_ops);


--
-- Name: assays_doc_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assays_doc_id ON assays USING btree (doc_id);


--
-- Name: assays_relationship_type; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assays_relationship_type ON assays USING btree (relationship_type);


--
-- Name: assays_relationship_type_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assays_relationship_type_like ON assays USING btree (relationship_type varchar_pattern_ops);


--
-- Name: assays_src_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assays_src_id ON assays USING btree (src_id);


--
-- Name: assays_tid; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX assays_tid ON assays USING btree (tid);


--
-- Name: atc_classification_level5_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX atc_classification_level5_like ON atc_classification USING btree (level5 varchar_pattern_ops);


--
-- Name: binding_sites_tid; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX binding_sites_tid ON binding_sites USING btree (tid);


--
-- Name: biotherapeutic_components_component_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX biotherapeutic_components_component_id ON biotherapeutic_components USING btree (component_id);


--
-- Name: biotherapeutic_components_molregno; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX biotherapeutic_components_molregno ON biotherapeutic_components USING btree (molregno);


--
-- Name: chembl_id_lookup_chembl_id_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX chembl_id_lookup_chembl_id_like ON chembl_id_lookup USING btree (chembl_id varchar_pattern_ops);


--
-- Name: component_class_component_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX component_class_component_id ON component_class USING btree (component_id);


--
-- Name: component_class_protein_class_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX component_class_protein_class_id ON component_class USING btree (protein_class_id);


--
-- Name: component_domains_component_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX component_domains_component_id ON component_domains USING btree (component_id);


--
-- Name: component_domains_domain_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX component_domains_domain_id ON component_domains USING btree (domain_id);


--
-- Name: component_sequences_accession_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX component_sequences_accession_like ON component_sequences USING btree (accession varchar_pattern_ops);


--
-- Name: component_synonyms_component_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX component_synonyms_component_id ON component_synonyms USING btree (component_id);


--
-- Name: compound_properties_alogp; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_properties_alogp ON compound_properties USING btree (alogp);


--
-- Name: compound_properties_hba; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_properties_hba ON compound_properties USING btree (hba);


--
-- Name: compound_properties_hbd; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_properties_hbd ON compound_properties USING btree (hbd);


--
-- Name: compound_properties_mw_freebase; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_properties_mw_freebase ON compound_properties USING btree (mw_freebase);


--
-- Name: compound_properties_num_ro5_violations; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_properties_num_ro5_violations ON compound_properties USING btree (num_ro5_violations);


--
-- Name: compound_properties_psa; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_properties_psa ON compound_properties USING btree (psa);


--
-- Name: compound_properties_rtb; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_properties_rtb ON compound_properties USING btree (rtb);


--
-- Name: compound_records_compound_key; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_records_compound_key ON compound_records USING btree (compound_key);


--
-- Name: compound_records_compound_key_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_records_compound_key_like ON compound_records USING btree (compound_key varchar_pattern_ops);


--
-- Name: compound_records_doc_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_records_doc_id ON compound_records USING btree (doc_id);


--
-- Name: compound_records_molregno; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_records_molregno ON compound_records USING btree (molregno);


--
-- Name: compound_records_src_compound_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_records_src_compound_id ON compound_records USING btree (src_compound_id);


--
-- Name: compound_records_src_compound_id_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_records_src_compound_id_like ON compound_records USING btree (src_compound_id varchar_pattern_ops);


--
-- Name: compound_records_src_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_records_src_id ON compound_records USING btree (src_id);


--
-- Name: compound_structures_standard_inchi_key; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_structures_standard_inchi_key ON compound_structures USING btree (standard_inchi_key);


--
-- Name: compound_structures_standard_inchi_key_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX compound_structures_standard_inchi_key_like ON compound_structures USING btree (standard_inchi_key varchar_pattern_ops);


--
-- Name: curation_lookup_curated_by_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX curation_lookup_curated_by_like ON curation_lookup USING btree (curated_by varchar_pattern_ops);


--
-- Name: data_validity_lookup_data_validity_comment_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX data_validity_lookup_data_validity_comment_like ON data_validity_lookup USING btree (data_validity_comment varchar_pattern_ops);


--
-- Name: defined_daily_dose_atc_code; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX defined_daily_dose_atc_code ON defined_daily_dose USING btree (atc_code);


--
-- Name: defined_daily_dose_atc_code_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX defined_daily_dose_atc_code_like ON defined_daily_dose USING btree (atc_code varchar_pattern_ops);


--
-- Name: django_session_expire_date; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX django_session_expire_date ON django_session USING btree (expire_date);


--
-- Name: django_session_session_key_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX django_session_session_key_like ON django_session USING btree (session_key varchar_pattern_ops);


--
-- Name: docs_chembl_id_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX docs_chembl_id_like ON docs USING btree (chembl_id varchar_pattern_ops);


--
-- Name: docs_issue; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX docs_issue ON docs USING btree (issue);


--
-- Name: docs_issue_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX docs_issue_like ON docs USING btree (issue varchar_pattern_ops);


--
-- Name: docs_journal; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX docs_journal ON docs USING btree (journal);


--
-- Name: docs_journal_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX docs_journal_like ON docs USING btree (journal varchar_pattern_ops);


--
-- Name: docs_volume; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX docs_volume ON docs USING btree (volume);


--
-- Name: docs_volume_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX docs_volume_like ON docs USING btree (volume varchar_pattern_ops);


--
-- Name: docs_year; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX docs_year ON docs USING btree (year);


--
-- Name: drug_mechanism_action_type; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX drug_mechanism_action_type ON drug_mechanism USING btree (action_type);


--
-- Name: drug_mechanism_action_type_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX drug_mechanism_action_type_like ON drug_mechanism USING btree (action_type varchar_pattern_ops);


--
-- Name: drug_mechanism_molregno; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX drug_mechanism_molregno ON drug_mechanism USING btree (molregno);


--
-- Name: drug_mechanism_record_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX drug_mechanism_record_id ON drug_mechanism USING btree (record_id);


--
-- Name: drug_mechanism_site_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX drug_mechanism_site_id ON drug_mechanism USING btree (site_id);


--
-- Name: drug_mechanism_tid; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX drug_mechanism_tid ON drug_mechanism USING btree (tid);


--
-- Name: formulations_molregno; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX formulations_molregno ON formulations USING btree (molregno);


--
-- Name: formulations_product_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX formulations_product_id ON formulations USING btree (product_id);


--
-- Name: formulations_product_id_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX formulations_product_id_like ON formulations USING btree (product_id varchar_pattern_ops);


--
-- Name: formulations_record_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX formulations_record_id ON formulations USING btree (record_id);


--
-- Name: fps2_atombv_idx; Type: INDEX; Schema: public; Owner: vagrant; Tablespace: 
--

CREATE INDEX fps2_atombv_idx ON fps2_rdkit USING gist (atombv);


--
-- Name: fps2_ffp2_idx; Type: INDEX; Schema: public; Owner: vagrant; Tablespace: 
--

CREATE INDEX fps2_ffp2_idx ON fps2_rdkit USING gist (ffp2);


--
-- Name: fps2_layfp_idx; Type: INDEX; Schema: public; Owner: vagrant; Tablespace: 
--

CREATE INDEX fps2_layfp_idx ON fps2_rdkit USING gist (layeredfp);


--
-- Name: fps2_maccsfp_idx; Type: INDEX; Schema: public; Owner: vagrant; Tablespace: 
--

CREATE INDEX fps2_maccsfp_idx ON fps2_rdkit USING gist (maccsfp);


--
-- Name: fps2_mfp2_idx; Type: INDEX; Schema: public; Owner: vagrant; Tablespace: 
--

CREATE INDEX fps2_mfp2_idx ON fps2_rdkit USING gist (mfp2);


--
-- Name: fps2_rdkfp_idx; Type: INDEX; Schema: public; Owner: vagrant; Tablespace: 
--

CREATE INDEX fps2_rdkfp_idx ON fps2_rdkit USING gist (rdkfp);


--
-- Name: fps2_ttbv_idx; Type: INDEX; Schema: public; Owner: vagrant; Tablespace: 
--

CREATE INDEX fps2_ttbv_idx ON fps2_rdkit USING gist (torsionbv);


--
-- Name: fps_ffp2_idx; Type: INDEX; Schema: public; Owner: vagrant; Tablespace: 
--

CREATE INDEX fps_ffp2_idx ON fps_rdkit USING gist (ffp2);


--
-- Name: fps_mfp2_idx; Type: INDEX; Schema: public; Owner: vagrant; Tablespace: 
--

CREATE INDEX fps_mfp2_idx ON fps_rdkit USING gist (mfp2);


--
-- Name: fps_ttbv_idx; Type: INDEX; Schema: public; Owner: vagrant; Tablespace: 
--

CREATE INDEX fps_ttbv_idx ON fps_rdkit USING gist (torsionbv);


--
-- Name: mechanism_refs_mec_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX mechanism_refs_mec_id ON mechanism_refs USING btree (mec_id);


--
-- Name: molecule_atc_classification_level5; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX molecule_atc_classification_level5 ON molecule_atc_classification USING btree (level5);


--
-- Name: molecule_atc_classification_level5_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX molecule_atc_classification_level5_like ON molecule_atc_classification USING btree (level5 varchar_pattern_ops);


--
-- Name: molecule_atc_classification_molregno; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX molecule_atc_classification_molregno ON molecule_atc_classification USING btree (molregno);


--
-- Name: molecule_dictionary_chembl_id_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX molecule_dictionary_chembl_id_like ON molecule_dictionary USING btree (chembl_id varchar_pattern_ops);


--
-- Name: molecule_dictionary_max_phase; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX molecule_dictionary_max_phase ON molecule_dictionary USING btree (max_phase);


--
-- Name: molecule_dictionary_pref_name; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX molecule_dictionary_pref_name ON molecule_dictionary USING btree (pref_name);


--
-- Name: molecule_dictionary_pref_name_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX molecule_dictionary_pref_name_like ON molecule_dictionary USING btree (pref_name varchar_pattern_ops);


--
-- Name: molecule_dictionary_therapeutic_flag; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX molecule_dictionary_therapeutic_flag ON molecule_dictionary USING btree (therapeutic_flag);


--
-- Name: molecule_hierarchy_active_molregno; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX molecule_hierarchy_active_molregno ON molecule_hierarchy USING btree (active_molregno);


--
-- Name: molecule_hierarchy_parent_molregno; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX molecule_hierarchy_parent_molregno ON molecule_hierarchy USING btree (parent_molregno);


--
-- Name: molecule_synonyms_molregno; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX molecule_synonyms_molregno ON molecule_synonyms USING btree (molregno);


--
-- Name: molecule_synonyms_res_stem_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX molecule_synonyms_res_stem_id ON molecule_synonyms USING btree (res_stem_id);


--
-- Name: parameter_type_parameter_type_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX parameter_type_parameter_type_like ON parameter_type USING btree (parameter_type varchar_pattern_ops);


--
-- Name: predicted_binding_domains_activity_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX predicted_binding_domains_activity_id ON predicted_binding_domains USING btree (activity_id);


--
-- Name: predicted_binding_domains_site_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX predicted_binding_domains_site_id ON predicted_binding_domains USING btree (site_id);


--
-- Name: products_product_id_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX products_product_id_like ON products USING btree (product_id varchar_pattern_ops);


--
-- Name: protein_class_synonyms_protein_class_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX protein_class_synonyms_protein_class_id ON protein_class_synonyms USING btree (protein_class_id);


--
-- Name: protein_family_classification_protein_class_desc_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX protein_family_classification_protein_class_desc_like ON protein_family_classification USING btree (protein_class_desc varchar_pattern_ops);


--
-- Name: rdkit_mol_idx; Type: INDEX; Schema: public; Owner: vagrant; Tablespace: 
--

CREATE INDEX rdkit_mol_idx ON compound_mols USING gist (ctab);


--
-- Name: relationship_type_relationship_type_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX relationship_type_relationship_type_like ON relationship_type USING btree (relationship_type varchar_pattern_ops);


--
-- Name: research_companies_res_stem_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX research_companies_res_stem_id ON research_companies USING btree (res_stem_id);


--
-- Name: research_stem_research_stem_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX research_stem_research_stem_like ON research_stem USING btree (research_stem varchar_pattern_ops);


--
-- Name: site_components_component_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX site_components_component_id ON site_components USING btree (component_id);


--
-- Name: site_components_domain_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX site_components_domain_id ON site_components USING btree (domain_id);


--
-- Name: site_components_site_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX site_components_site_id ON site_components USING btree (site_id);


--
-- Name: target_components_component_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX target_components_component_id ON target_components USING btree (component_id);


--
-- Name: target_components_tid; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX target_components_tid ON target_components USING btree (tid);


--
-- Name: target_dictionary_chembl_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX target_dictionary_chembl_id ON target_dictionary USING btree (chembl_id);


--
-- Name: target_dictionary_chembl_id_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX target_dictionary_chembl_id_like ON target_dictionary USING btree (chembl_id varchar_pattern_ops);


--
-- Name: target_dictionary_organism; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX target_dictionary_organism ON target_dictionary USING btree (organism);


--
-- Name: target_dictionary_organism_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX target_dictionary_organism_like ON target_dictionary USING btree (organism varchar_pattern_ops);


--
-- Name: target_dictionary_pref_name; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX target_dictionary_pref_name ON target_dictionary USING btree (pref_name);


--
-- Name: target_dictionary_pref_name_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX target_dictionary_pref_name_like ON target_dictionary USING btree (pref_name varchar_pattern_ops);


--
-- Name: target_dictionary_target_type; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX target_dictionary_target_type ON target_dictionary USING btree (target_type);


--
-- Name: target_dictionary_target_type_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX target_dictionary_target_type_like ON target_dictionary USING btree (target_type varchar_pattern_ops);


--
-- Name: target_dictionary_tax_id; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX target_dictionary_tax_id ON target_dictionary USING btree (tax_id);


--
-- Name: target_relations_related_tid; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX target_relations_related_tid ON target_relations USING btree (related_tid);


--
-- Name: target_relations_tid; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX target_relations_tid ON target_relations USING btree (tid);


--
-- Name: target_type_target_type_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX target_type_target_type_like ON target_type USING btree (target_type varchar_pattern_ops);


--
-- Name: version_name_like; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX version_name_like ON version USING btree (name varchar_pattern_ops);


--
-- Name: ws_cache_expires; Type: INDEX; Schema: public; Owner: chembl; Tablespace: 
--

CREATE INDEX ws_cache_expires ON ws_cache USING btree (expires);


--
-- Name: activities_assay_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY activities
    ADD CONSTRAINT activities_assay_id_fkey FOREIGN KEY (assay_id) REFERENCES assays(assay_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: activities_data_validity_comment_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY activities
    ADD CONSTRAINT activities_data_validity_comment_fkey FOREIGN KEY (data_validity_comment) REFERENCES data_validity_lookup(data_validity_comment) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: activities_doc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY activities
    ADD CONSTRAINT activities_doc_id_fkey FOREIGN KEY (doc_id) REFERENCES docs(doc_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: activities_molregno_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY activities
    ADD CONSTRAINT activities_molregno_fkey FOREIGN KEY (molregno) REFERENCES molecule_dictionary(molregno) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: activities_record_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY activities
    ADD CONSTRAINT activities_record_id_fkey FOREIGN KEY (record_id) REFERENCES compound_records(record_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: assay_parameters_assay_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY assay_parameters
    ADD CONSTRAINT assay_parameters_assay_id_fkey FOREIGN KEY (assay_id) REFERENCES assays(assay_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: assay_parameters_parameter_type_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY assay_parameters
    ADD CONSTRAINT assay_parameters_parameter_type_fkey FOREIGN KEY (parameter_type) REFERENCES parameter_type(parameter_type) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: assays_assay_type_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY assays
    ADD CONSTRAINT assays_assay_type_fkey FOREIGN KEY (assay_type) REFERENCES assay_type(assay_type) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: assays_cell_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY assays
    ADD CONSTRAINT assays_cell_id_fkey FOREIGN KEY (cell_id) REFERENCES cell_dictionary(cell_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: assays_chembl_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY assays
    ADD CONSTRAINT assays_chembl_id_fkey FOREIGN KEY (chembl_id) REFERENCES chembl_id_lookup(chembl_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: assays_confidence_score_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY assays
    ADD CONSTRAINT assays_confidence_score_fkey FOREIGN KEY (confidence_score) REFERENCES confidence_score_lookup(confidence_score) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: assays_curated_by_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY assays
    ADD CONSTRAINT assays_curated_by_fkey FOREIGN KEY (curated_by) REFERENCES curation_lookup(curated_by) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: assays_doc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY assays
    ADD CONSTRAINT assays_doc_id_fkey FOREIGN KEY (doc_id) REFERENCES docs(doc_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: assays_relationship_type_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY assays
    ADD CONSTRAINT assays_relationship_type_fkey FOREIGN KEY (relationship_type) REFERENCES relationship_type(relationship_type) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: assays_src_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY assays
    ADD CONSTRAINT assays_src_id_fkey FOREIGN KEY (src_id) REFERENCES source(src_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: assays_tid_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY assays
    ADD CONSTRAINT assays_tid_fkey FOREIGN KEY (tid) REFERENCES target_dictionary(tid) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: auth_group_permissions_permission_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY auth_group_permissions
    ADD CONSTRAINT auth_group_permissions_permission_id_fkey FOREIGN KEY (permission_id) REFERENCES auth_permission(id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: auth_permission_content_type_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY auth_permission
    ADD CONSTRAINT auth_permission_content_type_id_fkey FOREIGN KEY (content_type_id) REFERENCES django_content_type(id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: auth_user_groups_group_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY auth_user_groups
    ADD CONSTRAINT auth_user_groups_group_id_fkey FOREIGN KEY (group_id) REFERENCES auth_group(id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: auth_user_user_permissions_permission_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY auth_user_user_permissions
    ADD CONSTRAINT auth_user_user_permissions_permission_id_fkey FOREIGN KEY (permission_id) REFERENCES auth_permission(id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: binding_sites_tid_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY binding_sites
    ADD CONSTRAINT binding_sites_tid_fkey FOREIGN KEY (tid) REFERENCES target_dictionary(tid) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: biotherapeutic_components_component_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY biotherapeutic_components
    ADD CONSTRAINT biotherapeutic_components_component_id_fkey FOREIGN KEY (component_id) REFERENCES bio_component_sequences(component_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: biotherapeutic_components_molregno_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY biotherapeutic_components
    ADD CONSTRAINT biotherapeutic_components_molregno_fkey FOREIGN KEY (molregno) REFERENCES biotherapeutics(molregno) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: biotherapeutics_molregno_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY biotherapeutics
    ADD CONSTRAINT biotherapeutics_molregno_fkey FOREIGN KEY (molregno) REFERENCES molecule_dictionary(molregno) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: component_class_component_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY component_class
    ADD CONSTRAINT component_class_component_id_fkey FOREIGN KEY (component_id) REFERENCES component_sequences(component_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: component_class_protein_class_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY component_class
    ADD CONSTRAINT component_class_protein_class_id_fkey FOREIGN KEY (protein_class_id) REFERENCES protein_classification(protein_class_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: component_domains_component_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY component_domains
    ADD CONSTRAINT component_domains_component_id_fkey FOREIGN KEY (component_id) REFERENCES component_sequences(component_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: component_domains_domain_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY component_domains
    ADD CONSTRAINT component_domains_domain_id_fkey FOREIGN KEY (domain_id) REFERENCES domains(domain_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: component_synonyms_component_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY component_synonyms
    ADD CONSTRAINT component_synonyms_component_id_fkey FOREIGN KEY (component_id) REFERENCES component_sequences(component_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: compound_properties_molregno_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY compound_properties
    ADD CONSTRAINT compound_properties_molregno_fkey FOREIGN KEY (molregno) REFERENCES molecule_dictionary(molregno) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: compound_records_doc_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY compound_records
    ADD CONSTRAINT compound_records_doc_id_fkey FOREIGN KEY (doc_id) REFERENCES docs(doc_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: compound_records_molregno_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY compound_records
    ADD CONSTRAINT compound_records_molregno_fkey FOREIGN KEY (molregno) REFERENCES molecule_dictionary(molregno) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: compound_records_src_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY compound_records
    ADD CONSTRAINT compound_records_src_id_fkey FOREIGN KEY (src_id) REFERENCES source(src_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: compound_structures_molregno_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY compound_structures
    ADD CONSTRAINT compound_structures_molregno_fkey FOREIGN KEY (molregno) REFERENCES molecule_dictionary(molregno) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: defined_daily_dose_atc_code_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY defined_daily_dose
    ADD CONSTRAINT defined_daily_dose_atc_code_fkey FOREIGN KEY (atc_code) REFERENCES atc_classification(level5) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: docs_chembl_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY docs
    ADD CONSTRAINT docs_chembl_id_fkey FOREIGN KEY (chembl_id) REFERENCES chembl_id_lookup(chembl_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: drug_mechanism_action_type_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY drug_mechanism
    ADD CONSTRAINT drug_mechanism_action_type_fkey FOREIGN KEY (action_type) REFERENCES action_type(action_type) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: drug_mechanism_molregno_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY drug_mechanism
    ADD CONSTRAINT drug_mechanism_molregno_fkey FOREIGN KEY (molregno) REFERENCES molecule_dictionary(molregno) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: drug_mechanism_record_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY drug_mechanism
    ADD CONSTRAINT drug_mechanism_record_id_fkey FOREIGN KEY (record_id) REFERENCES compound_records(record_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: drug_mechanism_site_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY drug_mechanism
    ADD CONSTRAINT drug_mechanism_site_id_fkey FOREIGN KEY (site_id) REFERENCES binding_sites(site_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: drug_mechanism_tid_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY drug_mechanism
    ADD CONSTRAINT drug_mechanism_tid_fkey FOREIGN KEY (tid) REFERENCES target_dictionary(tid) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: formulations_molregno_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY formulations
    ADD CONSTRAINT formulations_molregno_fkey FOREIGN KEY (molregno) REFERENCES molecule_dictionary(molregno) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: formulations_product_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY formulations
    ADD CONSTRAINT formulations_product_id_fkey FOREIGN KEY (product_id) REFERENCES products(product_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: formulations_record_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY formulations
    ADD CONSTRAINT formulations_record_id_fkey FOREIGN KEY (record_id) REFERENCES compound_records(record_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: group_id_refs_id_f4b32aac; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY auth_group_permissions
    ADD CONSTRAINT group_id_refs_id_f4b32aac FOREIGN KEY (group_id) REFERENCES auth_group(id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: ligand_eff_activity_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY ligand_eff
    ADD CONSTRAINT ligand_eff_activity_id_fkey FOREIGN KEY (activity_id) REFERENCES activities(activity_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: mechanism_refs_mec_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY mechanism_refs
    ADD CONSTRAINT mechanism_refs_mec_id_fkey FOREIGN KEY (mec_id) REFERENCES drug_mechanism(mec_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: molecule_atc_classification_level5_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY molecule_atc_classification
    ADD CONSTRAINT molecule_atc_classification_level5_fkey FOREIGN KEY (level5) REFERENCES atc_classification(level5) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: molecule_atc_classification_molregno_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY molecule_atc_classification
    ADD CONSTRAINT molecule_atc_classification_molregno_fkey FOREIGN KEY (molregno) REFERENCES molecule_dictionary(molregno) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: molecule_dictionary_chembl_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY molecule_dictionary
    ADD CONSTRAINT molecule_dictionary_chembl_id_fkey FOREIGN KEY (chembl_id) REFERENCES chembl_id_lookup(chembl_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: molecule_hierarchy_active_molregno_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY molecule_hierarchy
    ADD CONSTRAINT molecule_hierarchy_active_molregno_fkey FOREIGN KEY (active_molregno) REFERENCES molecule_dictionary(molregno) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: molecule_hierarchy_molregno_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY molecule_hierarchy
    ADD CONSTRAINT molecule_hierarchy_molregno_fkey FOREIGN KEY (molregno) REFERENCES molecule_dictionary(molregno) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: molecule_hierarchy_parent_molregno_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY molecule_hierarchy
    ADD CONSTRAINT molecule_hierarchy_parent_molregno_fkey FOREIGN KEY (parent_molregno) REFERENCES molecule_dictionary(molregno) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: molecule_synonyms_molregno_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY molecule_synonyms
    ADD CONSTRAINT molecule_synonyms_molregno_fkey FOREIGN KEY (molregno) REFERENCES molecule_dictionary(molregno) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: molecule_synonyms_res_stem_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY molecule_synonyms
    ADD CONSTRAINT molecule_synonyms_res_stem_id_fkey FOREIGN KEY (res_stem_id) REFERENCES research_stem(res_stem_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: predicted_binding_domains_activity_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY predicted_binding_domains
    ADD CONSTRAINT predicted_binding_domains_activity_id_fkey FOREIGN KEY (activity_id) REFERENCES activities(activity_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: predicted_binding_domains_site_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY predicted_binding_domains
    ADD CONSTRAINT predicted_binding_domains_site_id_fkey FOREIGN KEY (site_id) REFERENCES binding_sites(site_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: protein_class_synonyms_protein_class_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY protein_class_synonyms
    ADD CONSTRAINT protein_class_synonyms_protein_class_id_fkey FOREIGN KEY (protein_class_id) REFERENCES protein_classification(protein_class_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: research_companies_res_stem_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY research_companies
    ADD CONSTRAINT research_companies_res_stem_id_fkey FOREIGN KEY (res_stem_id) REFERENCES research_stem(res_stem_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: site_components_component_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY site_components
    ADD CONSTRAINT site_components_component_id_fkey FOREIGN KEY (component_id) REFERENCES component_sequences(component_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: site_components_domain_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY site_components
    ADD CONSTRAINT site_components_domain_id_fkey FOREIGN KEY (domain_id) REFERENCES domains(domain_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: site_components_site_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY site_components
    ADD CONSTRAINT site_components_site_id_fkey FOREIGN KEY (site_id) REFERENCES binding_sites(site_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: target_components_component_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY target_components
    ADD CONSTRAINT target_components_component_id_fkey FOREIGN KEY (component_id) REFERENCES component_sequences(component_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: target_components_tid_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY target_components
    ADD CONSTRAINT target_components_tid_fkey FOREIGN KEY (tid) REFERENCES target_dictionary(tid) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: target_dictionary_chembl_id_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY target_dictionary
    ADD CONSTRAINT target_dictionary_chembl_id_fkey FOREIGN KEY (chembl_id) REFERENCES chembl_id_lookup(chembl_id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: target_dictionary_target_type_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY target_dictionary
    ADD CONSTRAINT target_dictionary_target_type_fkey FOREIGN KEY (target_type) REFERENCES target_type(target_type) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: target_relations_related_tid_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY target_relations
    ADD CONSTRAINT target_relations_related_tid_fkey FOREIGN KEY (related_tid) REFERENCES target_dictionary(tid) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: target_relations_tid_fkey; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY target_relations
    ADD CONSTRAINT target_relations_tid_fkey FOREIGN KEY (tid) REFERENCES target_dictionary(tid) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: user_id_refs_id_40c41112; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY auth_user_groups
    ADD CONSTRAINT user_id_refs_id_40c41112 FOREIGN KEY (user_id) REFERENCES auth_user(id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: user_id_refs_id_4dc23c39; Type: FK CONSTRAINT; Schema: public; Owner: chembl
--

ALTER TABLE ONLY auth_user_user_permissions
    ADD CONSTRAINT user_id_refs_id_4dc23c39 FOREIGN KEY (user_id) REFERENCES auth_user(id) DEFERRABLE INITIALLY DEFERRED;


--
-- Name: public; Type: ACL; Schema: -; Owner: postgres
--

REVOKE ALL ON SCHEMA public FROM PUBLIC;
REVOKE ALL ON SCHEMA public FROM postgres;
GRANT ALL ON SCHEMA public TO postgres;
GRANT ALL ON SCHEMA public TO PUBLIC;


--
-- Name: action_type; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE action_type FROM PUBLIC;
REVOKE ALL ON TABLE action_type FROM cbh_chembl_usr;
GRANT ALL ON TABLE action_type TO cbh_chembl_usr;


--
-- Name: activities; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE activities FROM PUBLIC;
REVOKE ALL ON TABLE activities FROM cbh_chembl_usr;
GRANT ALL ON TABLE activities TO cbh_chembl_usr;


--
-- Name: activity_stds_lookup; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE activity_stds_lookup FROM PUBLIC;
REVOKE ALL ON TABLE activity_stds_lookup FROM cbh_chembl_usr;
GRANT ALL ON TABLE activity_stds_lookup TO cbh_chembl_usr;


--
-- Name: assay_parameters; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE assay_parameters FROM PUBLIC;
REVOKE ALL ON TABLE assay_parameters FROM cbh_chembl_usr;
GRANT ALL ON TABLE assay_parameters TO cbh_chembl_usr;


--
-- Name: assay_type; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE assay_type FROM PUBLIC;
REVOKE ALL ON TABLE assay_type FROM cbh_chembl_usr;
GRANT ALL ON TABLE assay_type TO cbh_chembl_usr;


--
-- Name: assays; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE assays FROM PUBLIC;
REVOKE ALL ON TABLE assays FROM cbh_chembl_usr;
GRANT ALL ON TABLE assays TO cbh_chembl_usr;


--
-- Name: atc_classification; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE atc_classification FROM PUBLIC;
REVOKE ALL ON TABLE atc_classification FROM cbh_chembl_usr;
GRANT ALL ON TABLE atc_classification TO cbh_chembl_usr;


--
-- Name: binding_sites; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE binding_sites FROM PUBLIC;
REVOKE ALL ON TABLE binding_sites FROM cbh_chembl_usr;
GRANT ALL ON TABLE binding_sites TO cbh_chembl_usr;


--
-- Name: bio_component_sequences; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE bio_component_sequences FROM PUBLIC;
REVOKE ALL ON TABLE bio_component_sequences FROM cbh_chembl_usr;
GRANT ALL ON TABLE bio_component_sequences TO cbh_chembl_usr;


--
-- Name: biotherapeutic_components; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE biotherapeutic_components FROM PUBLIC;
REVOKE ALL ON TABLE biotherapeutic_components FROM cbh_chembl_usr;
GRANT ALL ON TABLE biotherapeutic_components TO cbh_chembl_usr;


--
-- Name: biotherapeutics; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE biotherapeutics FROM PUBLIC;
REVOKE ALL ON TABLE biotherapeutics FROM cbh_chembl_usr;
GRANT ALL ON TABLE biotherapeutics TO cbh_chembl_usr;


--
-- Name: cell_dictionary; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE cell_dictionary FROM PUBLIC;
REVOKE ALL ON TABLE cell_dictionary FROM cbh_chembl_usr;
GRANT ALL ON TABLE cell_dictionary TO cbh_chembl_usr;


--
-- Name: chembl_id_lookup; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE chembl_id_lookup FROM PUBLIC;
REVOKE ALL ON TABLE chembl_id_lookup FROM cbh_chembl_usr;
GRANT ALL ON TABLE chembl_id_lookup TO cbh_chembl_usr;


--
-- Name: component_class; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE component_class FROM PUBLIC;
REVOKE ALL ON TABLE component_class FROM cbh_chembl_usr;
GRANT ALL ON TABLE component_class TO cbh_chembl_usr;


--
-- Name: component_domains; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE component_domains FROM PUBLIC;
REVOKE ALL ON TABLE component_domains FROM cbh_chembl_usr;
GRANT ALL ON TABLE component_domains TO cbh_chembl_usr;


--
-- Name: component_sequences; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE component_sequences FROM PUBLIC;
REVOKE ALL ON TABLE component_sequences FROM cbh_chembl_usr;
GRANT ALL ON TABLE component_sequences TO cbh_chembl_usr;


--
-- Name: component_synonyms; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE component_synonyms FROM PUBLIC;
REVOKE ALL ON TABLE component_synonyms FROM cbh_chembl_usr;
GRANT ALL ON TABLE component_synonyms TO cbh_chembl_usr;


--
-- Name: compound_properties; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE compound_properties FROM PUBLIC;
REVOKE ALL ON TABLE compound_properties FROM cbh_chembl_usr;
GRANT ALL ON TABLE compound_properties TO cbh_chembl_usr;


--
-- Name: compound_records; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE compound_records FROM PUBLIC;
REVOKE ALL ON TABLE compound_records FROM cbh_chembl_usr;
GRANT ALL ON TABLE compound_records TO cbh_chembl_usr;


--
-- Name: compound_structures; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE compound_structures FROM PUBLIC;
REVOKE ALL ON TABLE compound_structures FROM cbh_chembl_usr;
GRANT ALL ON TABLE compound_structures TO cbh_chembl_usr;


--
-- Name: confidence_score_lookup; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE confidence_score_lookup FROM PUBLIC;
REVOKE ALL ON TABLE confidence_score_lookup FROM cbh_chembl_usr;
GRANT ALL ON TABLE confidence_score_lookup TO cbh_chembl_usr;


--
-- Name: curation_lookup; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE curation_lookup FROM PUBLIC;
REVOKE ALL ON TABLE curation_lookup FROM cbh_chembl_usr;
GRANT ALL ON TABLE curation_lookup TO cbh_chembl_usr;


--
-- Name: data_validity_lookup; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE data_validity_lookup FROM PUBLIC;
REVOKE ALL ON TABLE data_validity_lookup FROM cbh_chembl_usr;
GRANT ALL ON TABLE data_validity_lookup TO cbh_chembl_usr;


--
-- Name: defined_daily_dose; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE defined_daily_dose FROM PUBLIC;
REVOKE ALL ON TABLE defined_daily_dose FROM cbh_chembl_usr;
GRANT ALL ON TABLE defined_daily_dose TO cbh_chembl_usr;


--
-- Name: docs; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE docs FROM PUBLIC;
REVOKE ALL ON TABLE docs FROM cbh_chembl_usr;
GRANT ALL ON TABLE docs TO cbh_chembl_usr;


--
-- Name: domains; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE domains FROM PUBLIC;
REVOKE ALL ON TABLE domains FROM cbh_chembl_usr;
GRANT ALL ON TABLE domains TO cbh_chembl_usr;


--
-- Name: drug_mechanism; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE drug_mechanism FROM PUBLIC;
REVOKE ALL ON TABLE drug_mechanism FROM cbh_chembl_usr;
GRANT ALL ON TABLE drug_mechanism TO cbh_chembl_usr;


--
-- Name: formulations; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE formulations FROM PUBLIC;
REVOKE ALL ON TABLE formulations FROM cbh_chembl_usr;
GRANT ALL ON TABLE formulations TO cbh_chembl_usr;


--
-- Name: ligand_eff; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE ligand_eff FROM PUBLIC;
REVOKE ALL ON TABLE ligand_eff FROM cbh_chembl_usr;
GRANT ALL ON TABLE ligand_eff TO cbh_chembl_usr;


--
-- Name: mechanism_refs; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE mechanism_refs FROM PUBLIC;
REVOKE ALL ON TABLE mechanism_refs FROM cbh_chembl_usr;
GRANT ALL ON TABLE mechanism_refs TO cbh_chembl_usr;


--
-- Name: molecule_atc_classification; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE molecule_atc_classification FROM PUBLIC;
REVOKE ALL ON TABLE molecule_atc_classification FROM cbh_chembl_usr;
GRANT ALL ON TABLE molecule_atc_classification TO cbh_chembl_usr;


--
-- Name: molecule_dictionary; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE molecule_dictionary FROM PUBLIC;
REVOKE ALL ON TABLE molecule_dictionary FROM cbh_chembl_usr;
GRANT ALL ON TABLE molecule_dictionary TO cbh_chembl_usr;


--
-- Name: molecule_hierarchy; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE molecule_hierarchy FROM PUBLIC;
REVOKE ALL ON TABLE molecule_hierarchy FROM cbh_chembl_usr;
GRANT ALL ON TABLE molecule_hierarchy TO cbh_chembl_usr;


--
-- Name: molecule_synonyms; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE molecule_synonyms FROM PUBLIC;
REVOKE ALL ON TABLE molecule_synonyms FROM cbh_chembl_usr;
GRANT ALL ON TABLE molecule_synonyms TO cbh_chembl_usr;


--
-- Name: octmp_summary; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE octmp_summary FROM PUBLIC;
REVOKE ALL ON TABLE octmp_summary FROM cbh_chembl_usr;
GRANT ALL ON TABLE octmp_summary TO cbh_chembl_usr;


--
-- Name: organism_class; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE organism_class FROM PUBLIC;
REVOKE ALL ON TABLE organism_class FROM cbh_chembl_usr;
GRANT ALL ON TABLE organism_class TO cbh_chembl_usr;


--
-- Name: parameter_type; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE parameter_type FROM PUBLIC;
REVOKE ALL ON TABLE parameter_type FROM cbh_chembl_usr;
GRANT ALL ON TABLE parameter_type TO cbh_chembl_usr;


--
-- Name: predicted_binding_domains; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE predicted_binding_domains FROM PUBLIC;
REVOKE ALL ON TABLE predicted_binding_domains FROM cbh_chembl_usr;
GRANT ALL ON TABLE predicted_binding_domains TO cbh_chembl_usr;


--
-- Name: products; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE products FROM PUBLIC;
REVOKE ALL ON TABLE products FROM cbh_chembl_usr;
GRANT ALL ON TABLE products TO cbh_chembl_usr;


--
-- Name: protein_class_synonyms; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE protein_class_synonyms FROM PUBLIC;
REVOKE ALL ON TABLE protein_class_synonyms FROM cbh_chembl_usr;
GRANT ALL ON TABLE protein_class_synonyms TO cbh_chembl_usr;


--
-- Name: protein_classification; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE protein_classification FROM PUBLIC;
REVOKE ALL ON TABLE protein_classification FROM cbh_chembl_usr;
GRANT ALL ON TABLE protein_classification TO cbh_chembl_usr;


--
-- Name: protein_family_classification; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE protein_family_classification FROM PUBLIC;
REVOKE ALL ON TABLE protein_family_classification FROM cbh_chembl_usr;
GRANT ALL ON TABLE protein_family_classification TO cbh_chembl_usr;


--
-- Name: relationship_type; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE relationship_type FROM PUBLIC;
REVOKE ALL ON TABLE relationship_type FROM cbh_chembl_usr;
GRANT ALL ON TABLE relationship_type TO cbh_chembl_usr;


--
-- Name: research_companies; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE research_companies FROM PUBLIC;
REVOKE ALL ON TABLE research_companies FROM cbh_chembl_usr;
GRANT ALL ON TABLE research_companies TO cbh_chembl_usr;


--
-- Name: research_stem; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE research_stem FROM PUBLIC;
REVOKE ALL ON TABLE research_stem FROM cbh_chembl_usr;
GRANT ALL ON TABLE research_stem TO cbh_chembl_usr;


--
-- Name: site_components; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE site_components FROM PUBLIC;
REVOKE ALL ON TABLE site_components FROM cbh_chembl_usr;
GRANT ALL ON TABLE site_components TO cbh_chembl_usr;


--
-- Name: source; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE source FROM PUBLIC;
REVOKE ALL ON TABLE source FROM cbh_chembl_usr;
GRANT ALL ON TABLE source TO cbh_chembl_usr;


--
-- Name: target_components; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE target_components FROM PUBLIC;
REVOKE ALL ON TABLE target_components FROM cbh_chembl_usr;
GRANT ALL ON TABLE target_components TO cbh_chembl_usr;


--
-- Name: target_dictionary; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE target_dictionary FROM PUBLIC;
REVOKE ALL ON TABLE target_dictionary FROM cbh_chembl_usr;
GRANT ALL ON TABLE target_dictionary TO cbh_chembl_usr;


--
-- Name: target_relations; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE target_relations FROM PUBLIC;
REVOKE ALL ON TABLE target_relations FROM cbh_chembl_usr;
GRANT ALL ON TABLE target_relations TO cbh_chembl_usr;


--
-- Name: target_type; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE target_type FROM PUBLIC;
REVOKE ALL ON TABLE target_type FROM cbh_chembl_usr;
GRANT ALL ON TABLE target_type TO cbh_chembl_usr;


--
-- Name: usan_stems; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE usan_stems FROM PUBLIC;
REVOKE ALL ON TABLE usan_stems FROM cbh_chembl_usr;
GRANT ALL ON TABLE usan_stems TO cbh_chembl_usr;


--
-- Name: version; Type: ACL; Schema: public; Owner: chembl
--

REVOKE ALL ON TABLE version FROM PUBLIC;
REVOKE ALL ON TABLE version FROM cbh_chembl_usr;
GRANT ALL ON TABLE version TO cbh_chembl_usr;


--
-- PostgreSQL database dump complete
--

