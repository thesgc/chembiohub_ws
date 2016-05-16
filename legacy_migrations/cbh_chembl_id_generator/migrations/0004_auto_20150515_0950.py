# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    '''Migration to add the functions necessary in order to generate ids at the database level'''

    dependencies = [
        ('cbh_chembl_id_generator',
         '0003_remove_cbhcompoundid_original_project_key'),
    ]

    operations = [
        migrations.RunSQL("""
create or replace function random_int_string(length integer) returns text as 
$$
declare
  ---
  --- Set the chars to numbers
  ---
  chars text[] := '{0,1,2,3,4,5,6,7,8,9}';
  result text := '';
  i integer := 0;
begin
  if length < 0 then
    raise exception 'Given length cannot be less than 0';
  end if;
  for i in 1..length loop
    result := result || chars[1+random()*(array_length(chars, 1)-1)];
  end loop;
  return result;
end;
$$ language plpgsql;



create or replace function random_string(length integer) returns text as 
$$
declare
  ---
  --- Set the chars to capital letters
  ---
  chars text[] := '{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z}';
  result text := '';
  i integer := 0;
begin
  if length < 0 then
    raise exception 'Given length cannot be less than 0';
  end if;
  for i in 1..length loop
    result := result || chars[1+random()*(array_length(chars, 1)-1)];
  end loop;
  return result;
end;
$$ language plpgsql;


CREATE or replace FUNCTION make_new_id(structure_keyv char, id_key char, original_installation_keyv char)  
 RETURNS SETOF cbh_chembl_id_generator_cbhcompoundid
  VOLATILE
AS $$
---
--- Declare a rowtype so all input variables can be used as are
--- Expects the structure_keyv to be empty if this is a new secret compound
---
declare r cbh_chembl_id_generator_cbhcompoundid%rowtype;
declare new_id varchar(12):= id_key || random_string(3) || random_int_string(2) || random_string(2);
BEGIN
    LOOP
    for r in
---
---Empty values are not updated so it will create a new value. If the compound is
---Secret and you want a new batch you must send the old ID as the structure key - 
---secret compounds are stored with the UOX id in both columns
---
            UPDATE cbh_chembl_id_generator_cbhcompoundid SET current_batch_id=current_batch_id+1 WHERE structure_key = structure_keyv and structure_key != '' returning * 
            loop return next r;  
    END LOOP;
    IF found THEN
        RETURN;
    END IF;
    BEGIN
    IF structure_keyv='' THEN
        structure_keyv := new_id;
    END IF;
    RETURN QUERY INSERT INTO cbh_chembl_id_generator_cbhcompoundid ("structure_key", "assigned_id","current_batch_id",  "original_installation_key" ) values 
      (structure_keyv, new_id , 1, original_installation_keyv) returning *;
    RETURN;
    EXCEPTION WHEN unique_violation THEN
    END;
    END LOOP;
END;
$$ 
LANGUAGE plpgsql;
            """)
    ]
