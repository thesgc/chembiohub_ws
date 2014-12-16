'use strict';

/**
 * @ngdoc overview
 * @name ngChemApp
 * @description
 * # ngChemApp
 *
 * Main module of the application.
 */
angular.module('ngChemApp', [
    'ngAnimate',
    'ngAria',
    'ngCookies',
    'ngMessages',
    'ngResource',
    'ngGrid',
    'ngSanitize',
    'ngTouch',
    'ui.router',
    'ui.bootstrap'
  ])
  .config(function ($stateProvider, $urlRouterProvider) {
      
          // catch all route
    // send users to the form page 
    $urlRouterProvider.otherwise('/demo/intro');


      $stateProvider
          
        // HOME STATES AND NESTED VIEWS ========================================
        .state('demo', {
            url: '/demo',
            templateUrl: 'views/start.html',
            controller: 'DemoCtrl'
        })
        
        // nested states 
        // each of these sections will have their own view
        // url will be nested (/form/profile)
        .state('demo.intro', {
            url: '/intro',
            templateUrl: 'views/demo-intro.html',
            controller: function($rootScope) {
              $rootScope.step = 0;
              $rootScope.dynamic = 0;
              applyTicks("intro");
            }
        })

        .state('demo.add', {
            url: '/add',
            templateUrl: 'views/demo-add.html',
            controller: function($rootScope) {
              $rootScope.step = 1;
              $rootScope.dynamic = 41;
              applyTicks("add");
            }
        })

        .state('demo.map', {
            url: '/map',
            templateUrl: 'views/demo-map.html',
            controller: function($rootScope) {
              $rootScope.step = 2;
              $rootScope.dynamic = 58;
              applyTicks("map");
            }
        })
        
        // url will be /form/interests
        .state('demo.validate', {
            url: '/validate',
            templateUrl: 'views/demo-validate.html',
            controller: function($rootScope) {
              $rootScope.step = 3;
              $rootScope.dynamic = 74;
              applyTicks("validate");
            }
        })
        
        // url will be /form/payment
        .state('demo.finish', {
            url: '/finish',
            templateUrl: 'views/demo-finish.html',
            controller: function($rootScope) {
              $rootScope.step = 4;
              $rootScope.dynamic = 90.5;
              applyTicks("finish");
            }
        });
        

  });
