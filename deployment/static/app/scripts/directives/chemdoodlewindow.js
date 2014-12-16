'use strict';

/**
 * @ngdoc directive
 * @name ngChemApp.directive:chemdoodleWindow
 * @description
 * # chemdoodleWindow
 */
angular.module('ngChemApp')
  .directive('chemdoodleWindow', function () {
    return {
      template: '<div class="col-xs-12" id="chemdoodle-holder"><div>{{localMolFile}}</div><canvas id="chemdoodle" ng-click="saveMol()"></canvas></div>',
      restrict: 'E',
      scope:{'sketchMolfile':'=sketchMolfile' },
      link: function postLink(scope, element, attrs) {
      	
        //scope watching for a change to the scope value
        //jquery watching for a mouseup event to trigger the change in scope value.
        element.bind('mousedown', function(){
        	scope.localMolFile = ChemDoodle.writeMOL(element.getMolecule());
        });
        var cd_width = jQuery('#chemdoodle-holder').width();
        element = new ChemDoodle.SketcherCanvas('chemdoodle', cd_width, 300, {oneMolecule:true});
        element.repaint();
        
      },
      controller: function($scope) {
      	$scope.localMolFile = "";
      	$scope.saveMol = function() {
      		$scope.sketchMolfile = btoa("\n\n" + $scope.localMolFile.split("ichemlabs.com")[1]);
      	}
      }
    };
  });