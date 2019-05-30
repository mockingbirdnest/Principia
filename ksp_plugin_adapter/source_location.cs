// Pretend that we have the caller file path and line number attributes, even
// though they do not exist in .NET 3.5; since this is actually substituted
// at compile time it will work as intended.

namespace System.Runtime.CompilerServices {

[AttributeUsage(AttributeTargets.Parameter, Inherited = false)]
public sealed class CallerFilePathAttribute : Attribute {}

[AttributeUsage(AttributeTargets.Parameter, Inherited = false)]
public sealed class CallerLineNumberAttribute : Attribute {}

}  // namespace System.Runtime.CompilerServices
