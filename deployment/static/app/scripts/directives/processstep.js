'use strict';

/**
 * @ngdoc directive
 * @name ngChemApp.directive:processStep
 * @description
 * # processStep
 */
angular.module('ngChemApp')
  .directive('processStep', function () {
    return {
      template: '',
      restrict: 'E',
      link: function postLink(scope, element, attrs) {
        scope.getImageURL();
      },
      controller: ['$scope', function($scope) {
            $scope.getImageURL = function() {
                switch ($scope.status ){
                    case "no_change_req" :
                    $scope.imageURL = "images/skip_next2.svg";
                    break;
                    case "error" :
                    $scope.imageURL = "images/error2.svg";
                    break;
                    case "processing" :
                    $scope.imageURL = "images/loader.gif";
                    break;
                    case "changed" :
                    $scope.imageURL = "images/warning2.svg";
                    break;

                };

            }
    }],
      scope: {
        startSmiles : "=",
        endSmiles : "=",
        status : "=",
        name : "="
      }
    };
  });
