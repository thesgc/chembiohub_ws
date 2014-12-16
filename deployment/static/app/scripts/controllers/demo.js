'use strict';

/**
 * @ngdoc function
 * @name ngChemApp.controller:DemoCtrl
 * @description
 * # DemoCtrl
 * Controller of the ngChemApp
 */
var app = angular.module('ngChemApp');

app.controller('DemoCtrl', function ($scope, $rootScope, $state, ChEMBLFactory) {

	$scope.formData = {};

    $scope.sketchMolfile = "";

	$scope.myData = [ ];

    $rootScope.step = 0;
    $rootScope.totalSteps = 4;
    $rootScope.dynamic = 0;
    $rootScope.max = 100;

	$scope.alerts = [
        /*{ type: 'danger', msg: 'Oh snap! Change a few things up and try submitting again.' }*/
	];

	$scope.parsed_input = [ {"smiles": "COc1ccc2[C@@H]3[C@H](COc2c1)C(C)(C)OC4=C3C(=O)C(=O)C5=C4OC(C)(C)[C@H]6COc7cc(OC)ccc7[C@@H]56", "chemblId": "CHEMBL446858", "passesRuleOfThree": "No", "molecularWeight": 544.59, "molecularFormula": "C32H32O8", "acdLogp": 7.67, "knownDrug": "No", "stdInChiKey": "GHBOEFUAGSHXPO-UWXQAFAOSA-N", "species": null, "synonyms": null, "medChemFriendly": "Yes", "rotatableBonds": 2, "acdBasicPka": null, "alogp": 3.63, "preferredCompoundName": null, "numRo5Violations": 1, "acdLogd": 7.67, "acdAcidicPka": null},
                      {"smiles": "CN1C(=O)N(C)c2ncn(C)c2C1=O", "chemblId": "CHEMBL113", "passesRuleOfThree": "Yes", "molecularWeight": 194.19, "molecularFormula": "C8H10N4O2", "acdLogp": -0.63, "knownDrug": "Yes", "stdInChiKey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N", "species": "NEUTRAL", "synonyms": "Coffeine,SID104171124,SID11110908,Methyltheobromine,Caffeine,Cafcit,NoDoz Caplets and Chewable Tablets,SID124879556,Theine", "medChemFriendly": "Yes", "rotatableBonds": 0, "acdBasicPka": 0.52, "alogp": -0.10, "preferredCompoundName": "CAFFEINE", "numRo5Violations": 0, "acdLogd": -0.63, "acdAcidicPka": null} ];
    $scope.parsed_input.map(function(d){d.smiles = btoa(d.smiles)});
	$scope.input_string = "";

	$scope.gridOptions = { data: 'parsed_input',
                           showColumnMenu:true,
                            enableColumnReordering:true,
                            columnDefs: [{ field: "smiles", displayName: "Structure", cellTemplate:"img-template.html" },
                                        { field: "chemblId", displayName: "Chembl ID" },
                                        { field: "molecularWeight", displayName: "Mol Weight" },
                                        { field: "knownDrug", displayName: "Known Drug" },
                                        { field: "stdInChiKey", displayName: "Std InChi Key" },
                                        { field: "acdLogp", displayName: "acdLogp" }]
                          };

	$scope.parseInput = function (){
		//take string
		//split on delims
		//try and pattern match
		//work out what sort of input it is.
		//for now push our test json to the Input

		//now try a delimited list of SMILES
		var splitted = splitInput(this.input_string);


		angular.each(splitted, function(idx,val) {
			
			//intervene here with Chembl calls using the ChEMBLFactory
        	//interrogate each val for its type, then call the appropriate ChEMBL service to retrieve data
        	//var returned_json = {}
        	if (val.type == 'InChi') {
        		//val.extraData = ChEMBLFactory.queryInChi();
        		$scope.parsed_input.push(new ChEMBLFactory(val.input_str, "inchi"));
        		
        	}
        	else if (val.type == 'InChi Key') {
        		//val.extraData = ChEMBLFactory.queryInChiKey();
        		$scope.parsed_input.push(new ChEMBLFactory(val.input_str, "inchikey"));
        		
        		
        	}
        	else if (val.type == 'SMILES') {
        		$scope.parsed_input = (new ChEMBLFactory(val.input_str, "smiles"));
        		
        		
        	}
        	else if (val.type == 'ChEMBL ID') {
        		$scope.parsed_input.push(new ChEMBLFactory(val.input_str, "chemblid"));
        		
        	}

        	//display errors
        	/*if (returned_json.resp == "danger") {
        		$scope.addAlert(returned_json.resp, returned_json.message );
        	}*/
			//$scope.parsed_input.push(val);	
			console.log($scope.parsed_input);
		});



		//$scope.parsed_input.push(splitted);
	};

    $scope.addAlert = function(message, alert_type){
      $scope.alerts.push({msg: message || 'There has been an error!', type : alert_type || "danger" });
    };

    $scope.closeAlert = function(index) {
      $scope.alerts.splice(index, 1);
    };


    $scope.validatedData = [
                                [
                                    {
                                        name: "Original Structure",
                                        smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                                        status: "start"
                                    },
                                    {
                                        name: "Metals and Salts",
                                        smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                                        status: "no_change_req"
                                    },
                                    {
                                        name: "Tautomers and Charges",
                                        smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                                        status: "no_change_req"
                                    },
                                    {
                                        name: "Isotopes and Stereochemistry",
                                        smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                                        status: "changed"
                                    },
                                    {
                                        name: "Normalised Structure",
                                        smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                                        status: "end"
                                    }
                                ],
                                [
                                    {
                                        name: "Original",
                                        smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                                        status: "start"
                                    },
                                    {
                                        name: "Metals and Salts",
                                        smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                                        status: "no_change_req"
                                    },
                                    {
                                        name: "Tautomers and Charges",
                                        smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                                        status: "no_change_req"
                                    },
                                    {
                                        name: "Isotopes and Stereochemistry",
                                        smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                                        status: "processing"
                                    },
                                    {
                                        name: "Complete",
                                        smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                                        status: "processing"
                                    }
                                ],
                                [
                                    {
                                        name: "Original",
                                        smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                                        status: "start"
                                    },
                                    {
                                        name: "Metals and Salts",
                                        smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                                        status: "no_change_req"
                                    },
                                    {
                                        name: "Tautomers and Charges",
                                        smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                                        status: "error"
                                    },
                                    {
                                        name: "Isotopes and Stereochemistry",
                                        smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                                        status: "error"
                                    },
                                    {
                                        name: "Complete",
                                        smiles: "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
                                        status: "error"
                                    }
                                ]
                            ];



});

app.controller('MessageCtrl', function ($scope){
	

});


