'use strict';

/**
 * @ngdoc function
 * @name ngChemApp.controller:MainCtrl
 * @description
 * # MainCtrl
 * Controller of the ngChemApp
 */
angular.module('ngChemApp')
  .controller('MainCtrl', function ($scope) {
    $scope.awesomeThings = [
      'HTML5 Boilerplate',
      'AngularJS',
      'Karma'
    ];
  });
