'use strict';

/**
 * @ngdoc function
 * @name ngChemApp.controller:AboutCtrl
 * @description
 * # AboutCtrl
 * Controller of the ngChemApp
 */
angular.module('ngChemApp')
  .controller('AboutCtrl', function ($scope) {
    $scope.awesomeThings = [
      'HTML5 Boilerplate',
      'AngularJS',
      'Karma'
    ];
  });
