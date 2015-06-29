

      google.load("visualization", "1.1", {packages:["wordtree",'corechart',"table"]});
      google.setOnLoadCallback(drawSimpleNodeChart);


      function drawSimpleNodeChart() {
        var nodeListData = new google.visualization.DataTable();
        nodeListData.addColumn('number', 'id');
        nodeListData.addColumn('string', 'childLabel');
        nodeListData.addColumn('number', 'parent');
        nodeListData.addColumn('number', 'size');
        nodeListData.addColumn('number', 'color');

        nodeListData.addRow([0, 'Life', -1, 1, 0]);

        nodeListData.addRow([1, 'Boreoeutheria', 0, 4, 0]);

        nodeListData.addRow([2, 'Euarchontoglires', 1, 3, 0]);
        nodeListData.addRow([3, 'Canis Familiaris', 1, 1, 2]);


        nodeListData.addRow([4, 'Homininae', 2, 2, 0]);
        nodeListData.addRow([5, 'Murinae', 2, 2, 0]);

        nodeListData.addRow([6, 'Mus Musculus', 5, 1, 2]);
        nodeListData.addRow([7, 'Ratus Norvegicus', 5, 1, 2]);

        nodeListData.addRow([8, 'Homo Sapiens', 4, 1, 2]);
        nodeListData.addRow([9, 'Pan Troglodytes', 4, 1, 2]);
        nodeListData.addRow([10, 'Gorilla Gorilla', 4, 1, 2]);



        var options = {
          colors: ['black', 'black', 'black'],
          wordtree: {
            maxFontSize: 3,
            format: 'explicit',
            fontName: 'Helvetica Neue',
            type: 'suffix',
            backgroundColor: '#a98'
          }
        };

        var wordtree = new google.visualization.WordTree(document.getElementById('wordtree_explicit'));
        wordtree.draw(nodeListData, options);

        google.visualization.events.addListener(wordtree, 'select', selectHandler);

        function selectHandler(e) {
            drawCharttaxon(wordtree.getSelection().word);
        }




      }

       function drawCharttaxon(taxon) {
       $.getJSON("lineage.json", function(data_taxon) {
        var data = new google.visualization.DataTable();
      data.addColumn('string', 'taxon');
      data.addColumn('number', 'FULL TAXON');
      data.addColumn('number', 'HP');
      data.addColumn('number', 'GHP');
      data.addColumn('number', 'MR');




    $.each(data_taxon['perfect_chart'], function (key, f) {
    if (key === taxon){
        data.addRows([
        [key,  data_taxon['perfect_chart'][key][0],data_taxon['perfect_chart'][key][1],data_taxon['perfect_chart'][key][2],data_taxon['perfect_chart'][key][3] ]

        ]);
        }

    })



        var options = {
          title: '% of perfect matches depending of taxon' + taxon,
          curveType: 'none',

        };

        var chart = new google.visualization.ColumnChart(document.getElementById('donut_taxon'));

        chart.draw(data, options);
        });
      }

