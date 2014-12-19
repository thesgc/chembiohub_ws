'use strict';

/**
 * @ngdoc directive
 * @name ngChemApp.directive:chemsvg
 * @description
 * # chemsvg
 */
angular.module('ngChemApp')
  .directive('chemsvg', function ($compile) {
    return {
      template: '<img ng-src="{{baseUrl}}{{encSmiles}}/{{size}}"></img>',
      restrict: 'E',
      link: function postLink(scope, element, attrs) {
        //element.html('this is the chemsvg directive');
        scope.getSmiles();
      },
      controller: ['$scope', function($scope) {
            $scope.getSmiles = function() {
                $scope.encSmiles = btoa($scope.smiles);
            };
        }],
      scope: {
        size:'=',
        smiles: '=',
        baseUrl: '=',
      }
    };
  });
