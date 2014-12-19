'use strict';

var app = angular.module('ngChemApp');

/** Factories for querying APIs **/
app.factory('ChEMBLFactory',[ '$http', function ($http) {
	
	//efine your web service location here.

	

	var ChEMBLFactory = function(input_str, input_type) {
		
		//initialize by storing the input string as a variable
		var self = this;
		//self.input_str = input_str;
		//self.input_type = input_type;
		var chembl_base_url = "http://dee.sgc.ox.ac.uk:8000/chemblws/compounds/";

		this.queryChembl = function(query_string) {
        	//send smiles off to ChEMBL web service - populate returned data into our table
        	var chembl_data = $http.get(query_string);
			chembl_data.then(function(response){
				//angular.extend(self, response.data);
				self.data = response.data.compounds[0];
			});
        };

		if(self.input_type == "inchi") {
			this.queryChembl(chembl_base_url + self.input_str + '.json');
		}
		else if(self.input_type == "inchikey") {
			this.queryChembl(chembl_base_url + 'stdinchikey/' + self.input_str + '.json');
		}
		else if(self.input_type == "smiles") {
			this.queryChembl(chembl_base_url + 'smiles/' + self.input_str + '.json');
		}
		else if(self.input_type == "chemblid") {
			this.queryChembl(chembl_base_url + self.input_str + '.json');
		}

	};

	return (ChEMBLFactory);
}]);