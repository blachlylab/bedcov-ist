/** bedcov via Interval Splay Trees

    mirror approach as in https://github.com/lh3/cgranges
*/

module bedcov_ist;

import std.algorithm : splitter;
import std.array : appender;
import std.conv : to;
import std.range : take;
import std.stdio;
import std.string : toStringz;

import core.stdc.stdio : printf;
import core.stdc.stdint;
import core.stdc.stdlib;    // atoi

import intervaltree;
version(avltree) import intervaltree.avltree;
version(splaytree) import intervaltree.splaytree;
version(iitree) import intervaltree.cgranges;

/// parse bed line (credit: lh3)
char *parse_bed(char *s, int32_t *st_, int32_t *en_)
{
	char * p, q, ctg = null;
	int32_t i, st = -1, en = -1;
	for (i = 0, p = q = s;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
			if (i == 0) ctg = p;
			else if (i == 1) st = atoi(p);
			else if (i == 2) en = atoi(p);
			++i, p = q + 1;
			if (c == '\0') break;
		}
	}
	*st_ = st, *en_ = en;
	return i >= 3? ctg : null;
}

int main(string[] args)
{
    version(iitree)
    {
        // ugly
        //IITree!BasicInterval treesByContig = IITree!(BasicInterval)(cr_init());

        auto tree = cr_init();
    }
    version(splaytree)
        IntervalSplayTree!(BasicInterval)*[string] treesByContig;
    version(avltree)
        IntervalAVLTree!(BasicInterval)*[string] treesByContig;

    if (args.length != 3) {
        writeln("Usage: bedcov-ist <loaded.bed> <streamed.bed>");
        return 1;
    }
    
    auto fi1 = File(args[1], "r");
    auto fi2 = File(args[2], "r");
    auto fields = appender!(char[][]);

    // Index the first BED file into an interval tree
    foreach(line; fi1.byLine())
    {
//        if (line.length > 0 && line[0] == '#') continue;
//        else if (line.length > 4 && (line[0..5] == "track" || line[0..5] == "brows")) continue;

        fields.clear();
        fields.put( line.splitter.take(3) );    // max BED3

        string contig = fields.data[0].idup;
        int start = fields.data[1].to!int;
        int end = fields.data[2].to!int;

        version(iitree)
        {
            //treesByContig.add(contig, BasicInterval(start, end));

            cr_add(tree, toStringz(contig), start, end, 0);
        }
        version(splaytree)
        {
            auto tree = treesByContig.require(contig, new IntervalSplayTree!BasicInterval);

            tree.insert(BasicInterval(start, end));
        }
        version(avltree)
        {
            auto tree = treesByContig.require(contig, new IntervalAVLTree!BasicInterval);

            uint cnt;
            auto node = new IntervalTreeNode!BasicInterval(BasicInterval(start,end));
            tree.insert(node, cnt);
        }
    }

    //version(iitree) treesByContig.index();
    version(iitree) cr_index(tree);

    // Now intersect the second bedfile
    version(iitree)
    {
        int64_t *b;
        int64_t m_b;
        int64_t n_b;
    }
    foreach(line; fi2.byLine())
    {
//        if (line.length > 0 && line[0] == '#') continue;
//        else if (line.length > 4 && (line[0..5] == "track" || line[0..5] == "brows")) continue;

        fields.clear();
        fields.put( line.splitter.take(3) );    // max BED3

        auto contig = fields.data[0];
        int start = fields.data[1].to!int;
        int end = fields.data[2].to!int;

        //contig[$] = '\0'; // range violation
        line[contig.length] = '\0';

        version(iitree)
        {
            //if (cr_get_ctg( treesByContig.cr, toStringz(contig) ) == -1) continue;
            if (cr_get_ctg( tree, contig.ptr ) == -1) continue;
        }
        else
        {
            auto tree = treesByContig.require(contig.idup, null);
            if (tree is null) continue;
        }

        version(iitree)
        {
            //auto o = treesByContig.findOverlapsWith(contig, BasicInterval(start, end));             // returns cr_intv_t*(s)
            n_b = cr_overlap(tree, contig.ptr, start, end, &b, &m_b);
            printf("%s\t%d\t%d\t%ld\n", contig.ptr, start, end, n_b);
        }
        else
        {
            auto o = treesByContig[contig].findOverlapsWith(BasicInterval(start, end));   // returns Node*(s)
            printf("%s\t%d\t%d\t%ld\n", contig.ptr, start, end, o.length);
        }


    }

    return 0;
}