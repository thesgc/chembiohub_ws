'use strict';

/**
 * @ngdoc service
 * @name ngChemApp.CBHCompoundBatch
 * @description
 * # CBHCompoundBatch
 * Factory in the ngChemApp.
 */
angular.module('ngChemApp')
  .factory('CBHCompoundBatch', ['$http', '$q' ,function ($http) {
    // Service logic
    // ...

    var urlBase = "/chemblws/cbh_compound_batches/";
    var CBHCompoundBatch = {};

    CBHCompoundBatch.validate = function(molfile) {

      var arr = window.location.href.split("/");
      var myUrl = arr[0] + "//" + arr[2] + urlBase;

      return $http.post( myUrl + "validate/", {ctab:molfile});
    };

    return CBHCompoundBatch;

  }]);
