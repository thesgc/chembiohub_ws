'use strict';

describe('Filter: titleCase', function () {

  // load the filter's module
  beforeEach(module('ngChemApp'));

  // initialize a new instance of the filter before each test
  var titleCase;
  beforeEach(inject(function ($filter) {
    titleCase = $filter('titleCase');
  }));

  it('should return the input prefixed with "titleCase filter:"', function () {
    var text = 'angularjs';
    expect(titleCase(text)).toBe('titleCase filter: ' + text);
  });

});
