'use strict';

describe('Directive: processStep', function () {

  // load the directive's module
  beforeEach(module('ngChemApp'));

  var element,
    scope;

  beforeEach(inject(function ($rootScope) {
    scope = $rootScope.$new();
  }));

  it('should make hidden element visible', inject(function ($compile) {
    element = angular.element('<process-step></process-step>');
    element = $compile(element)(scope);
    expect(element.text()).toBe('this is the processStep directive');
  }));
});
