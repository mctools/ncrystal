// Windows-only developer tool (not part of the NCrystal distribution). See
// Analyzer.csproj's top comment for what this is and why it exists.
//
// Usage: ncrystal_winprofile_analyzer <trace.etl> <process-name-substring> [pdb-search-dir]
//
// Prints a flat "self time by leaf function" histogram over every kernel
// CPU-sampling event belonging to the first matching process found in the
// trace, ordered by descending sample count. This is deliberately the
// simplest possible aggregation (leaf frame only, not a full
// inclusive-time call tree): enough to answer "which function is hot",
// with much less API surface to get wrong than a full call-tree
// reconstruction, given this cannot be tested against a real Windows ETW
// trace before being run for real in CI.
//
// Matches CPU-sample events by EventName=="PerfInfo/Sample" -- confirmed
// via a real trace's diagnostic dump (see the fallback further down) after
// two earlier guesses ("SampledProfile", then anything containing both
// "Sample" and "Profile") both matched nothing despite the trace clearly
// containing stack-walked samples (10938 PerfInfo/Sample events were
// found sitting right alongside Thread/CSwitch and Dispatcher/ReadyThread,
// which are unrelated scheduler events, not CPU samples). The fallback
// dump is kept in case a different WPR profile or OS version ever uses a
// different name again.

using Microsoft.Diagnostics.Tracing.Etlx;
using Microsoft.Diagnostics.Symbols;

if (args.Length < 2)
{
    Console.Error.WriteLine("usage: ncrystal_winprofile_analyzer <trace.etl> <process-name-substring> [pdb-search-dir]");
    return 1;
}
string etlFile = args[0];
string procNameFilter = args[1];
string? pdbDir = args.Length > 2 ? args[2] : null;

try
{
    Console.WriteLine($"Opening {etlFile} ...");
    using var traceLog = TraceLog.OpenOrConvert(etlFile);

    var proc = traceLog.Processes
        .Where(p => p.Name.Contains(procNameFilter, StringComparison.OrdinalIgnoreCase))
        .OrderByDescending(p => p.CPUMSec)
        .FirstOrDefault();

    if (proc == null)
    {
        Console.Error.WriteLine($"Could not find a process matching \"{procNameFilter}\" in the trace.");
        Console.WriteLine("Processes seen in trace (top 20 by CPU):");
        foreach (var p in traceLog.Processes.OrderByDescending(p => p.CPUMSec).Take(20))
            Console.WriteLine($"  {p.Name} (PID {p.ProcessID}) CPU={p.CPUMSec}ms");
        return 1;
    }

    Console.WriteLine($"Process: {proc.Name} (PID {proc.ProcessID}), total CPU={proc.CPUMSec}ms");

    var symPathStr = SymbolPath.MicrosoftSymbolServerPath;
    if (!string.IsNullOrEmpty(pdbDir))
        symPathStr = pdbDir + ";" + symPathStr;
    using var symReader = new SymbolReader(Console.Out, symPathStr);

    foreach (var mod in proc.LoadedModules)
    {
        try { traceLog.CodeAddresses.LookupSymbolsForModule(symReader, mod.ModuleFile); }
        catch (Exception ex) { Console.WriteLine($"  (symbol lookup failed for {mod.FilePath}: {ex.Message})"); }
    }

    var selfTime = new Dictionary<string, int>();
    var eventNameCountsWithStacks = new Dictionary<string, int>();
    int totalSamples = 0;
    foreach (var ev in traceLog.Events)
    {
        if (ev.ProcessID != proc.ProcessID)
            continue;
        var cs = ev.CallStack();
        if (cs == null)
            continue;
        //Track this regardless of event name, so that if the filter below
        //still ends up finding nothing, the printed breakdown makes the
        //actual event name to match on obvious from the log alone:
        var evName = ev.EventName ?? "?";
        eventNameCountsWithStacks[evName] = eventNameCountsWithStacks.TryGetValue(evName, out var ec) ? ec + 1 : 1;

        if (evName != "PerfInfo/Sample")
            continue;
        totalSamples++;
        var codeAddr = cs.CodeAddress;
        string name = codeAddr.FullMethodName;
        if (string.IsNullOrEmpty(name))
            name = $"0x{codeAddr.Address:x} in {codeAddr.ModuleName ?? "?"}";
        selfTime[name] = selfTime.TryGetValue(name, out var c) ? c + 1 : 1;
    }

    Console.WriteLine($"Total CPU samples for process: {totalSamples}");
    Console.WriteLine("Top 40 functions by self (leaf) sample count:");
    foreach (var kv in selfTime.OrderByDescending(kv => kv.Value).Take(40))
    {
        double pct = totalSamples > 0 ? 100.0 * kv.Value / totalSamples : 0;
        Console.WriteLine($"  {kv.Value,6} ({pct,5:F1}%)  {kv.Key}");
    }

    if (totalSamples == 0)
    {
        Console.WriteLine();
        Console.WriteLine("No \"PerfInfo/Sample\" events were found. Event names that DO have"
                          + " stacks attached for this process (so one of these is probably the"
                          + " real CPU-sample event to match on instead):");
        foreach (var kv in eventNameCountsWithStacks.OrderByDescending(kv => kv.Value))
            Console.WriteLine($"  {kv.Value,6}  {kv.Key}");
    }

    return 0;
}
catch (Exception ex)
{
    Console.Error.WriteLine("Analysis failed: " + ex);
    return 1;
}
