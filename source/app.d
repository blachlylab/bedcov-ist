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

import intervaltree;
import intervaltree.splaytree;
import intervaltree.cgranges;

int main(string[] args)
{
    version(iitree)
    {
        // ugly
        IITree!BasicInterval treesByContig = IITree!(BasicInterval)(cr_init());
    }
    else
    {
        IntervalSplayTree!(BasicInterval)*[string] treesByContig;
        IntervalSplayTree!(BasicInterval)*[string] treesByContig2;
    }

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
            treesByContig.add(contig, BasicInterval(start, end));
        }
        else
        {
            auto tree = treesByContig.require(contig, new IntervalSplayTree!BasicInterval);

            tree.insert(BasicInterval(start, end));
        }
    }

    version(iitree) treesByContig.index();

    // Now intersect the second bedfile
    foreach(line; fi2.byLine())
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
            if (cr_get_ctg( treesByContig.cr, toStringz(contig) ) == -1) continue;
        }
        else
        {
            auto tree = treesByContig.require(contig, null);
            if (tree is null) continue;
        }

        version(iitree)
            auto o = treesByContig.findOverlapsWith(contig, BasicInterval(start, end));             // returns cr_intv_t*(s)
        else
            auto o = treesByContig[contig].findOverlapsWith(BasicInterval(start, end));   // returns Node*(s)

        printf("%s\t%d\t%d\t%ld\n", toStringz(contig), start, end, o.length);

    }

    return 0;
}