angular.module('ngChemApp.directives',[])
.directive('draggable', function () {
	return {
	  restrict: 'A',
	  scope:false,
	//  templateUrl : 'templates/dimensions.html',
	  link: function postLink(scope, element, attrs) {

	      scope.$watch('metadata', function(metadata){
	      	if(!metadata.length) element.find('li').remove();
		      element.find('li').draggable({
		        connectToSortable:'.dimensions-container',
				helper : 'clone',
		        revert: 'invalid',
		        start : onStart
		      })
	     	})

		   	function onStart(e,ui){
		      ui.helper.width($(e.currentTarget).width())
		      ui.helper.css('z-index','100000');
		    }

	  }
	};
})
.directive('droppable', function(){
  return {
    restrict: 'A',
    scope:false,
    //scope: {} //defines input parameters as string literal(@) or javascript to be run later(=)
    link: function postLink(scope, element, attrs) {
    	element.sortable({
    		items : '> li',
	        connectWith: '.dimensions-container',
	        placeholder:'drop',
	        /*start: onStart,
	        update: onUpdate,
	        receive : onReceive,
	        remove: onRemove,
	        over: over,*/
	        tolerance:'intersect'
    	});

    	function over(e,ui){
	    	var dimension = ui.item.data().dimension,
	    	html = isValidType(dimension) ? '<i class="fa fa-arrow-circle-down breath-right"></i>Drop here' : '<i class="fa fa-times-circle breath-right"></i>Don\'t drop here'
	    	element.find('.drop').html(html);
        } 

		function onUpdate(e,ui){

			ui.item.find('.dimension-icon').remove();

	    	if (ui.item.find('span.remove').length == 0) {
	      	ui.item.append("<span class='remove pull-right'>&times;</span>")
	      }
	     	ui.item.find('span.remove').click(function(){  ui.item.remove(); onRemove(); });

	     	if (removeLast) {
	     		ui.item.remove();
	     		removeLast = false;
	     	}    	

	     	scope.value = values();
	     	scope.$apply();

	     	element.parent().css("overflow","hidden");

				var dimension = ui.item.data().dimension;
	     	ui.item.toggleClass("invalid", !isValidType(dimension))
	     	message();

	     	$rootScope.$broadcast("update");
	    }

	    //etc - functions to provide extra functionality and alter scope values attached to the droppable/sortable
    }
  };
});